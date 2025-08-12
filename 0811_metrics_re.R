# ===== 패키지 / 옵션 =====
library(sf); library(terra); library(landscapemetrics)
library(dplyr); library(purrr); library(tidyr)
sf_use_s2(FALSE); terraOptions(progress = 1, memfrac = 0.6)

root    <- "/Users/hyeseon/Downloads"
out_dir <- file.path(root, "outputs"); dir.create(out_dir, showWarnings = FALSE)
sur_dir <- file.path(out_dir, "mw_surfaces"); dir.create(sur_dir, showWarnings = FALSE)
crs5179 <- "EPSG:5179"

# ===== 1) 입력 로드 & 5179 통일 =====
aoi_sf   <- st_read(file.path(root, "boundary.gpkg"), quiet=TRUE) |> st_make_valid() |> st_transform(5179)
lulc_sf  <- st_read(file.path(root, "lulc_1m_whole.shp"), quiet=TRUE) |> st_make_valid() |> st_transform(5179)
stopifnot("L1_CODE" %in% names(lulc_sf))
aoi  <- aoi_sf
lulc <- vect(lulc_sf)

# ===== 2) LULC → 10 m 정수 래스터 =====
tmpl   <- rast(ext(vect(aoi)), resolution = 10, crs = crs5179)
lulc10 <- rasterize(lulc, tmpl, field = "L1_CODE", touches = TRUE, background = NA) |>
  crop(vect(aoi)) |> mask(vect(aoi)) |> as.int()
writeRaster(lulc10, file.path(out_dir, "lulc_10m.tif"),
            overwrite = TRUE, gdal = c("COMPRESS=LZW"))

# ===== 3) 100 m 격자 & 벡터 변환 =====
grid100 <- st_make_grid(aoi, cellsize = 100, what = "polygons") |>
  st_intersection(aoi) |> st_sf() |> mutate(grid_id = dplyr::row_number())
st_write(grid100, file.path(out_dir, "grid_100m.gpkg"), delete_dsn = TRUE)
grid_v <- vect(grid100)

# ===== 4) window_lsm 표면(300/500/700 m) 만들기 =====
win_m  <- c(300, 500, 700)
odd    <- function(n) ifelse(n %% 2 == 0, n + 1, n)
win_px <- odd(round(win_m / res(lulc10)[1])); names(win_px) <- paste0(win_m, "m")
mk_win <- function(n) matrix(1, nrow = n, ncol = n)
metric_map <- c(lsm_l_shdi = "SHDI", lsm_l_lpi = "LPI", lsm_l_contag = "CONTAG")

# ★ 핵심: window_lsm() 반환을 항상 SpatRaster로 정규화
as_spatr <- function(x){
  if (inherits(x, "SpatRaster")) return(x)
  if (inherits(x, "list")){
    x <- x[!vapply(x, is.null, logical(1))]
    if (length(x) == 1) return(rast(x[[1]]))
    return(do.call(c, lapply(x, rast)))  # 여러 레이어면 스택으로 결합
  }
  rast(x)
}

run_surface <- function(metric, wpx, tag){
  res <- window_lsm(
    landscape = lulc10,
    what      = metric,          # 하나씩만
    window    = mk_win(wpx),
    na.rm     = TRUE,
    progress  = TRUE
  )
  r <- as_spatr(res)             # ← list 방지
  nm <- paste0(metric_map[[metric]], "_", tag)
  names(r) <- nm
  f <- file.path(sur_dir, paste0(nm, ".tif"))
  writeRaster(r, f, overwrite = TRUE, gdal = c("COMPRESS=LZW"))
  rast(f)
}

surface_list <- purrr::imap(win_px, ~c(
  run_surface("lsm_l_shdi",   .x, .y),
  run_surface("lsm_l_lpi",    .x, .y),
  run_surface("lsm_l_contag", .x, .y)  # 클래스 1종 창에서는 CONTAG = NA(정상 경고)
))

# ==== 5) 표면 → 100 m 격자 평균 요약 (수정 버전) ====
lm_tbl <- purrr::imap(surface_list, ~{
  df <- terra::extract(.x, grid_v, fun = function(v) mean(v, na.rm = TRUE))
  tibble(grid_id = grid100$grid_id) |> bind_cols(df[ , -1])
}) |> purrr::reduce(left_join, by = "grid_id")


# (선택) CONTAG_* 의 NA를 100으로 치환(완전 응집으로 간주)
lm_tbl <- lm_tbl |> mutate(across(starts_with("CONTAG_"), ~ ifelse(is.na(.x), 100, .x)))

# ===== 6) 저장 =====
write.csv(lm_tbl, file.path(out_dir, "landscape_metrics_by_grid_100m.csv"), row.names = FALSE)
st_write(left_join(grid100, lm_tbl, by = "grid_id"),
         file.path(out_dir, "landscape_metrics_100m.gpkg"), delete_dsn = TRUE)

# (선택) 스칼로그램 요약
scalo <- tibble(
  scale = names(win_px),
  SHDI_mean   = sapply(names(win_px), \(s) mean(lm_tbl[[paste0("SHDI_",   s)]], na.rm = TRUE)),
  LPI_mean    = sapply(names(win_px), \(s) mean(lm_tbl[[paste0("LPI_",    s)]], na.rm = TRUE)),
  CONTAG_mean = sapply(names(win_px), \(s) mean(lm_tbl[[paste0("CONTAG_", s)]], na.rm = TRUE))
)
write.csv(scalo, file.path(out_dir, "scalogram_summary.csv"), row.names = FALSE)
