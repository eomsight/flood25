
######## 10m grain, 300-500-700m moving window !
# ===== 패키지 / 옵션 =====
library(sf); library(terra); library(landscapemetrics)
library(dplyr); library(purrr); library(tidyr)
sf_use_s2(FALSE); terraOptions(progress = 1, memfrac = 0.6)

root    <- "/Users/hyeseon/Downloads"
out_dir <- file.path(root, "outputs");      dir.create(out_dir,  showWarnings = FALSE, recursive = TRUE)
sur_dir <- file.path(out_dir, "mw_surfaces"); dir.create(sur_dir, showWarnings = FALSE, recursive = TRUE)
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
writeRaster(lulc10, file.path(out_dir, "lulc_10m.tif"), overwrite = TRUE, gdal = c("COMPRESS=LZW"))

# ===== 4) window_lsm 표면(300/500/700 m) 만들기 =====
win_m  <- c(300, 500, 700)
odd    <- function(n) ifelse(n %% 2 == 0, n + 1, n)
win_px <- odd(round(win_m / res(lulc10)[1])); names(win_px) <- paste0(win_m, "m")
mk_win <- function(n) matrix(1, nrow = n, ncol = n)

metric_map <- c(
  lsm_l_shdi   = "SHDI",
  lsm_l_lpi    = "LPI",
  lsm_l_contag = "CONTAG",
  lsm_l_ed     = "ED"            # NEW (ED)
)

as_spatr <- function(x){
  if (inherits(x, "SpatRaster")) return(x)
  if (inherits(x, "list")){
    x <- x[!vapply(x, is.null, logical(1))]
    if (length(x) == 1) return(rast(x[[1]]))
    return(do.call(c, lapply(x, rast)))
  }
  rast(x)
}

run_surface <- function(metric, wpx, tag){
  res <- window_lsm(
    landscape = lulc10,
    what      = metric,
    window    = mk_win(wpx),
    na.rm     = TRUE,
    progress  = TRUE
  )
  r <- as_spatr(res)
  nm <- paste0(metric_map[[metric]], "_", tag)
  names(r) <- nm
  f <- file.path(sur_dir, paste0(nm, ".tif"))
  writeRaster(r, f, overwrite = TRUE, gdal = c("COMPRESS=LZW"))
  rast(f)
}

surface_list <- purrr::imap(win_px, ~c(
  run_surface("lsm_l_shdi",   .x, .y),
  run_surface("lsm_l_lpi",    .x, .y),
  run_surface("lsm_l_contag", .x, .y),
  run_surface("lsm_l_ed",     .x, .y)      # NEW (ED)
))

        # 10 m 표면 한 스택(12밴드)
sr <- c(surface_list[["300m"]], surface_list[["500m"]], surface_list[["700m"]])
names(sr) <- c("SHDI_300m","LPI_300m","CONTAG_300m","ED_300m",     # NEW (ED)
               "SHDI_500m","LPI_500m","CONTAG_500m","ED_500m",     # NEW (ED)
               "SHDI_700m","LPI_700m","CONTAG_700m","ED_700m")     # NEW (ED)

# === 10 m 결과 저장(택1 또는 둘 다) ===
writeRaster(sr, file.path(out_dir, "landscape_metrics_10m.tif"), overwrite = TRUE)
writeRaster(sr, file.path(out_dir, "landscape_metrics_10m.gpkg"),
            filetype = "GPKG", overwrite = TRUE, gdal = c("RASTER_TABLE" = "metrics_10m"))



# sr (12밴드: SHDI/LPI/CONTAG/ED × 300/500/700m) 만든 직후에 아래 3줄 추가
dir <- out_dir                                         # 네가 아래에서 쓰는 dir를 out_dir로 맵핑
m   <- as.data.frame(sr)                               # 10 m 셀 단위로 각 밴드 값을 데이터프레임으로
m$grid_id <- seq_len(nrow(m)); m <- dplyr::relocate(m, grid_id)  # grid_id 추가






####### optimal scale 찾기
### Scalogram: 평균·분산 안정 구간(plateau) 찾기
# (패키지/로드/dir/m/g 읽는 부분은 네 코드 그대로)

# long 포맷
mm <- m |>
  tidyr::pivot_longer(-grid_id, names_to = "name", values_to = "val") |>
  tidyr::separate(name, into = c("metric","scale"), sep = "_(?=\\d)") |>
  dplyr::mutate(scale = as.integer(stringr::str_remove(scale, "m")))

# scalogram 통계
scalo <- mm %>%
  dplyr::group_by(metric, scale) %>%
  dplyr::summarise(mean = mean(val, na.rm=TRUE),
                   sd   = sd(val, na.rm=TRUE),
                   cv   = sd/mean,
                   .groups = "drop") %>%
  dplyr::arrange(metric, scale)

# 상대변화(%)
scalo_delta <- scalo %>%
  dplyr::group_by(metric) %>%
  dplyr::arrange(scale, .by_group=TRUE) %>%
  dplyr::mutate(rel_change = abs((mean - dplyr::lag(mean))/dplyr::lag(mean))*100) %>%
  dplyr::ungroup()

# (그림은 그대로)
ggplot2::ggplot(scalo, aes(scale, mean)) +
  ggplot2::geom_line() + ggplot2::geom_point() +
  ggplot2::facet_wrap(~metric, scales="free_y") +
  ggplot2::theme_bw() + ggplot2::labs(x="Observation scale (m)", y="Mean metric")

# plateau 판정 (그대로)
pick_plateau <- scalo_delta %>%
  dplyr::group_by(metric) %>%
  dplyr::summarise(
    d1 = rel_change[scale==500],  # 300→500
    d2 = rel_change[scale==700],  # 500→700
    plateau = (!is.na(d1) & !is.na(d2) & d1 < 2 & d2 < 2),
    pick_by_plateau = ifelse(plateau, 300L, NA_integer_),
    .groups = "drop"
  )

pick_plateau

### Redundancy: 스케일 간 상관 (그대로 쓰되 ED까지 호출)
cors <- function(pre){ stats::cor(m[grep(paste0("^",pre,"_"), names(m))],
                                  use="pairwise.complete.obs") }
cors("SHDI"); cors("LPI"); cors("CONTAG"); cors("ED")      # NEW

# ======= (NEW) 자동 최종 스케일 선택: plateau + 500~700 중복성 =======
get_cor57 <- function(metric){                              # NEW
  c500 <- paste0(metric,"_500m"); c700 <- paste0(metric,"_700m")
  if (!all(c(c500,c700) %in% names(m))) return(NA_real_)
  suppressWarnings(stats::cor(m[[c500]], m[[c700]], use="pairwise.complete.obs"))
}

rho_thr <- 0.95                                             # NEW
metrics <- unique(scalo$metric)                             # NEW
cor57_tbl <- dplyr::tibble(metric = metrics,
                           cor57  = sapply(metrics, get_cor57))  # NEW

pick_final <- pick_plateau %>%                              # NEW
  dplyr::left_join(cor57_tbl, by="metric") %>%
  dplyr::mutate(
    pick = dplyr::case_when(
      plateau ~ 300L,                                      # 1) plateau면 300
      !is.na(cor57) & cor57 >= rho_thr ~ 500L,             # 2) 500≈700면 500
      TRUE ~ ifelse(d2 < d1, 700L, 500L)                   # 3) 그외 d2<d1→700, 아니면 500
    )
  )

pick_final    # ← 최종 표: metric, d1, d2, plateau, cor57, pick

# ======= (선택) '선정된 스케일'만 추출해서 저장 =======           # NEW
sel_cols <- paste0(pick_final$metric, "_", pick_final$pick, "m")
readr::write_csv(m[, c("grid_id", sel_cols)],
                 file.path(dir, "landscape_metrics_selected_scale.csv"), na="")




####### 4metrics -> scale 500 값으로 다 통일해서 저장할게요
global_pick <- 500L
sel_cols_global <- paste0(c("SHDI","LPI","CONTAG","ED"), "_", global_pick, "m")
readr::write_csv(m[, c("grid_id", sel_cols_global)],
                 file.path(dir, paste0("landscape_metrics_global_", global_pick, "m.csv")),
                 na = "NA")




####### 경계라인 NA값만 버퍼 줘서 추가분석으로 빵꾸 메워보겠습니다
## ---- 0) 현재 500m 스택(기존 계산) ----
global_pick <- 500L
r_500 <- c(
  sr[[paste0("SHDI_",   global_pick, "m")]],
  sr[[paste0("LPI_",    global_pick, "m")]],
  sr[[paste0("CONTAG_", global_pick, "m")]],
  sr[[paste0("ED_",     global_pick, "m")]]
)
names(r_500) <- c("SHDI_500m","LPI_500m","CONTAG_500m","ED_500m")

## ---- 1) AOI를 500m 창 반경만큼 넉넉히 버퍼 → LULC 10m 재래스터화(빠름) ----
rad_m   <- ceiling(global_pick / 2)         # 500m 창의 반경 = 250m
aoi_buf <- sf::st_buffer(aoi_sf, rad_m)

tmpl_buf   <- terra::rast(terra::ext(terra::vect(aoi_buf)), resolution = 10, crs = crs5179)
lulc10_buf <- terra::rasterize(lulc, tmpl_buf, field = "L1_CODE", touches = TRUE, background = NA) |>
  terra::crop(terra::vect(aoi_buf)) |>
  terra::mask(terra::vect(aoi_buf)) |>
  terra::as.int()

## ---- 2) 500m만 이동창 재계산(4개 지표) ----
wpx500 <- odd(round(global_pick / terra::res(lulc10_buf)[1]))

run_surface_buf <- function(metric, wpx){
  res <- landscapemetrics::window_lsm(
    landscape = lulc10_buf,
    what      = metric,
    window    = mk_win(wpx),
    na.rm     = TRUE,
    progress  = TRUE
  )
  as_spatr(res)
}

re500 <- c(
  run_surface_buf("lsm_l_shdi",   wpx500),
  run_surface_buf("lsm_l_lpi",    wpx500),
  run_surface_buf("lsm_l_contag", wpx500),
  run_surface_buf("lsm_l_ed",     wpx500)
)
names(re500) <- names(r_500)

# 원래 AOI로 다시 잘라주기
re500 <- terra::crop(re500, terra::vect(aoi_sf)) |>
  terra::mask(terra::vect(aoi_sf))

## ---- 3) NA만 교체해서 '채워진' 최종 500m 만들기 ----
# cover: 첫 번째의 NA를 두 번째의 값으로 대체
r_500_filled <- terra::cover(r_500, re500)

## ---- 4) 저장(GeoTIFF + CSV) ----
terra::writeRaster(r_500_filled,
                   file.path(out_dir, "maps_500m_filled.tif"),
                   overwrite = TRUE, gdal = c("COMPRESS=LZW")
)

df500 <- as.data.frame(r_500_filled)
df500$grid_id <- seq_len(nrow(df500))
df500 <- dplyr::relocate(df500, grid_id)
readr::write_csv(df500,
                 file.path(out_dir, "landscape_metrics_global_500m_filled.csv"),
                 na = "NA"
)

## ---- (선택) 얼마나 채워졌는지 확인 ----
before_na <- colSums(is.na(as.data.frame(r_500)))
after_na  <- colSums(is.na(as.data.frame(r_500_filled)))
print(rbind(before_na, after_na))







######## CONTAG NA -> 100으로 처리할게요
library(readr); library(dplyr)

csv_in  <- file.path(out_dir, "landscape_metrics_global_500m_filled.csv")
csv_out <- file.path(out_dir, "landscape_metrics_global_500m_filled_re2.csv")  # 새 버전

df <- read_csv(csv_in, show_col_types = FALSE)

tol <- 1e-6
# “실질 단일-클래스” 판정: SHDI≈0 & ED≈0 & CONTAG=NA  (LPI는 보지 않음)
single_holes <- (is.finite(df$SHDI_500m) & df$SHDI_500m <= tol) &
  (is.finite(df$ED_500m)   & df$ED_500m   <= tol) &
  is.na(df$CONTAG_500m)

# 치환
df$CONTAG_500m[single_holes] <- 100

# 남은 NA 카운트 확인
message("CONTAG NA 남은 건수: ", sum(is.na(df$CONTAG_500m)))

write_csv(df, csv_out, na = "NA")



######## GeoTIFF - QGIS 용 출력
library(terra)
tif_in  <- file.path(out_dir, "maps_500m_filled.tif")
tif_out <- file.path(out_dir, "maps_500m_filled_re2.tif")

r <- rast(tif_in)
names(r) <- c("SHDI_500m","LPI_500m","CONTAG_500m","ED_500m")

tol <- 1e-6
single_holes <- (r[["SHDI_500m"]] <= tol) &
  (r[["ED_500m"]]   <= tol) &
  is.na(r[["CONTAG_500m"]])

r[["CONTAG_500m"]] <- ifel(single_holes, 100, r[["CONTAG_500m"]])

writeRaster(r, tif_out, overwrite=TRUE, gdal=c("COMPRESS=LZW"))



