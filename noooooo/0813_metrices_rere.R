# ===== 패키지 / 옵션 =====
library(sf); library(terra); library(landscapemetrics)
library(dplyr); library(purrr); library(tidyr)
sf_use_s2(FALSE); terraOptions(progress = 1, memfrac = 0.6)

root    <- "/Users/hyeseon/Downloads"
out_dir <- file.path(root, "outputs");      dir.create(out_dir,  showWarnings = FALSE, recursive = TRUE)   # << 변경
sur_dir <- file.path(out_dir, "mw_surfaces"); dir.create(sur_dir, showWarnings = FALSE, recursive = TRUE)  # << 변경
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

# ===== 3) (선택) 100 m 격자 생성은 남겨둬도 OK — 어차피 아래에서 안 씀 =====
# grid100 <- st_make_grid(aoi, cellsize = 100, what = "polygons") |>
#   st_intersection(aoi) |> st_sf() |> mutate(grid_id = dplyr::row_number())
# st_write(grid100, file.path(out_dir, "grid_100m.gpkg"), delete_dsn = TRUE)
# grid_v <- vect(grid100)

# ===== 4) window_lsm 표면(300/500/700 m) 만들기 =====
win_m  <- c(300, 500, 700)
odd    <- function(n) ifelse(n %% 2 == 0, n + 1, n)
win_px <- odd(round(win_m / res(lulc10)[1])); names(win_px) <- paste0(win_m, "m")
mk_win <- function(n) matrix(1, nrow = n, ncol = n)
metric_map <- c(lsm_l_shdi = "SHDI", lsm_l_lpi = "LPI", lsm_l_contag = "CONTAG")

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
  rast(f)  # 그대로 반환 (예전과 동일)
}

surface_list <- purrr::imap(win_px, ~c(
  run_surface("lsm_l_shdi",   .x, .y),
  run_surface("lsm_l_lpi",    .x, .y),
  run_surface("lsm_l_contag", .x, .y)
))



# 10 m 표면 9개를 한 스택으로 묶기
sr <- c(surface_list[["300m"]], surface_list[["500m"]], surface_list[["700m"]])
names(sr) <- c("SHDI_300m","LPI_300m","CONTAG_300m",
               "SHDI_500m","LPI_500m","CONTAG_500m",
               "SHDI_700m","LPI_700m","CONTAG_700m")

# === 10 m 결과 저장(택1 또는 둘 다) ===
# 1) GeoTIFF(권장, QGIS 호환 최고)
writeRaster(sr, file.path(out_dir, "landscape_metrics_10m.tif"), overwrite = TRUE)

# 2) GeoPackage(한 파일에 9밴드)
writeRaster(sr, file.path(out_dir, "landscape_metrics_10m.gpkg"),
            filetype = "GPKG", overwrite = TRUE,
            gdal = c("RASTER_TABLE" = "metrics_10m"))


list_lsm(level = "landscape")
print(list_lsm(level = "landscape"), n = Inf)

