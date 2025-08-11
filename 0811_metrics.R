install.packages("sf")
install.packages("terra")
install.packages("landscapemetrics")
install.packages("dplyr")
install.packages("purrr")

# ==== 패키지 ====
library(sf); library(terra); library(landscapemetrics)
library(dplyr); library(purrr)
sf_use_s2(FALSE); terraOptions(progress = 1)

# ==== 경로/설정 ====
root <- "/Users/hyeseon/Downloads"
crs5179 <- "EPSG:5179"
path_boundary <- file.path(root, "boundary.gpkg")
path_lulc_vec <- file.path(root, "lulc_1m_whole.shp")   # 필드: L1_CODE (7분류)
out_dir <- file.path(root, "outputs"); dir.create(out_dir, showWarnings = FALSE)

# ---- 헬퍼(5179로 표준화) ----
to_5179_sf <- function(p) st_transform(st_make_valid(st_read(p, quiet=TRUE)), 5179)
to_5179_vect <- function(p) vect(to_5179_sf(p))

# ==== 1) 경계/LULC 5179로 ====
aoi  <- to_5179_sf(path_boundary)
lulc <- to_5179_vect(path_lulc_vec)  # L1_CODE 사용
st_write(aoi, file.path(out_dir, "boundary_5179.gpkg"), delete_dsn=TRUE)

# ==== 2) LULC → 10 m 정수 래스터 ====
tmpl   <- rast(ext(vect(aoi)), resolution = 10, crs = crs5179)
lulc10 <- rasterize(lulc, tmpl, field = "L1_CODE", touches = TRUE, background = NA) |>
  crop(vect(aoi)) |> mask(vect(aoi))
writeRaster(lulc10, file.path(out_dir, "lulc_10m.tif"),
            overwrite=TRUE, gdal=c("COMPRESS=LZW"))

# ==== 3) 100 m 격자 ====
grid100 <- st_make_grid(aoi, cellsize=100, what="polygons") |>
  st_intersection(aoi) |> st_sf() |> mutate(grid_id = row_number())
st_write(grid100, file.path(out_dir, "grid_100m.gpkg"), delete_dsn=TRUE)

# ==== 4) 경관지표(10 m 그레인, 창 300/500/700 m) ====
# 창 길이 → 픽셀수(홀수)
win_m  <- c(300, 500, 700)
odd    <- function(n) ifelse(n %% 2 == 0, n+1, n)
win_px <- odd(round(win_m / res(lulc10)[1]))
names(win_px) <- paste0(win_m, "m")

# 창 픽셀 → 행렬
mk_win <- function(n) matrix(1, nrow = n, ncol = n)

calc_win <- function(w, tag){
  wmat <- mk_win(w)
  r <- window_lsm(
    as.int(lulc10),                                   # 범주형 보장
    what   = c("lsm_l_shdi","lsm_l_lpi","lsm_l_contag"),
    window = wmat,
    na.rm  = TRUE,
    progress = TRUE
  )
  # r에 실제로 존재하는 레이어 이름을 기반으로 suffix만 붙이기
  # (lsm_l_shdi → SHDI 등으로 깔끔히 바꾸고 싶으면 매핑 사용)
  nm_map <- c(lsm_l_shdi="SHDI", lsm_l_lpi="LPI", lsm_l_contag="CONTAG")
  curr   <- names(r)
  nice   <- ifelse(curr %in% names(nm_map), nm_map[curr], curr)
  names(r) <- paste0(nice, "_", tag)
  r
}

stacks <- purrr::imap(win_px, ~calc_win(.x, .y))


# ==== 5) 100 m 격자 요약 & 저장 ====
grid_v <- vect(grid100)
lm_tbl <- purrr::imap(stacks, ~{
  df <- extract(.x, grid_v, fun = mean, na.rm = TRUE)
  tibble(grid_id = grid100$grid_id) |> bind_cols(df[,-1])
}) |> purrr::reduce(left_join, by="grid_id")

write.csv(lm_tbl, file.path(out_dir, "landscape_metrics_by_grid_100m.csv"), row.names=FALSE)
st_write(left_join(grid100, lm_tbl, by="grid_id"),
         file.path(out_dir, "landscape_metrics_100m.gpkg"), delete_dsn=TRUE)

# (선택) 스칼로그램 요약
scalo <- map_dfr(names(win_px), \(s)
                 tibble(scale=s,
                        SHDI_mean   = mean(lm_tbl[[paste0("SHDI_", s)]],   na.rm=TRUE),
                        LPI_mean    = mean(lm_tbl[[paste0("LPI_",  s)]],   na.rm=TRUE),
                        CONTAG_mean = mean(lm_tbl[[paste0("CONTAG_", s)]], na.rm=TRUE))
)
write.csv(scalo, file.path(out_dir, "scalogram_summary.csv"), row.names=FALSE)

