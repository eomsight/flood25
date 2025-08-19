library(landscapemetrics)
list_lsm(level = "class")
print(list_lsm(level = "class"), n = Inf)

install.packages("future.apply")


# ==== 패키지 ====
library(sf); library(terra); library(landscapemetrics); library(readr)
sf_use_s2(FALSE); terraOptions(progress = 1, memfrac = 0.6)

# ==== 경로/좌표계 ====
root    <- "/Users/hyeseon/Downloads"
out_dir <- file.path(root, "outputs"); dir.create(out_dir, TRUE, TRUE)
crs5179 <- "EPSG:5179"

# ==== 입력 로드 ====
aoi_sf  <- st_read(file.path(root, "boundary.gpkg"), quiet=TRUE) |> st_make_valid() |> st_transform(5179)
lulc_sf <- st_read(file.path(root, "lulc_1m_whole.shp"), quiet=TRUE) |> st_make_valid() |> st_transform(5179)
stopifnot("L1_CODE" %in% names(lulc_sf))
lulc    <- terra::vect(lulc_sf)

# ==== AOI 255 m 버퍼 → 10 m 래스터화 (경계 보정) ====
rad_m    <- 255L                                    # 51px(=510 m) 창과 정합
aoi_buf  <- sf::st_buffer(aoi_sf, rad_m)
tmpl_buf <- terra::rast(terra::ext(terra::vect(aoi_buf)), resolution = 10, crs = crs5179)

lulc10_buf <- terra::rasterize(lulc, tmpl_buf, field = "L1_CODE",
                               touches = TRUE, background = NA) |>
  terra::crop(terra::vect(aoi_buf)) |>
  terra::mask(terra::vect(aoi_buf)) |>
  terra::as.int()

# ==== 원래 AOI 10 m 그리드의 '중심 셀들' (계산 대상) ====
tmpl10 <- terra::rast(terra::ext(terra::vect(aoi_sf)), resolution = 10, crs = crs5179)
mask10 <- terra::init(tmpl10, 1) |>
  terra::crop(terra::vect(aoi_sf)) |>
  terra::mask(terra::vect(aoi_sf))

cells <- which(!is.na(terra::values(mask10)))
xy    <- terra::xyFromCell(mask10, cells)

# ---- (옵션) 너무 느리면 stride를 2,5 등으로 키워 샘플링 가능 ----
stride <- 1L                      # 1=모든 10m 셀(권장), 2=20m 간격, ...
cells <- cells[seq(1, length(cells), by = stride)]
xy    <- xy[seq(1, nrow(xy),     by = stride), , drop = FALSE]

# ==== 출력 벡터 준비 ====
n <- length(cells)
ed_vec  <- rep(NA_real_, n)
lpi_vec <- rep(NA_real_, n)

# ==== 한 지점 계산 함수 (클래스=100, ED & LPI) ====
calc_one <- function(x, y){
  ext_win <- terra::ext(x - rad_m, x + rad_m, y - rad_m, y + rad_m)  # 51px(510 m)
  tile <- terra::crop(lulc10_buf, ext_win, snap = "out")
  if (all(is.na(terra::values(tile)))) return(c(NA_real_, NA_real_))
  m <- suppressWarnings(
    landscapemetrics::calculate_lsm(
      landscape = tile,
      level     = "class",
      what      = c("lsm_c_ed","lsm_c_lpi"),
      classes   = 100,
      directions= 8,
      progress  = FALSE
    )
  )
  ed  <- m$value[m$metric == "lsm_c_ed"];  ed  <- if (length(ed)==0)  NA_real_ else ed
  lpi <- m$value[m$metric == "lsm_c_lpi"]; lpi <- if (length(lpi)==0) NA_real_ else lpi
  c(ed, lpi)   # units: ED=m/ha, LPI=%
}

# ==== 메인 루프(순차; 안정 우선) ====
pb <- txtProgressBar(min = 0, max = n, style = 3)
for (i in seq_len(n)) {
  vals <- calc_one(xy[i,1], xy[i,2])
  ed_vec[i]  <- vals[1]
  lpi_vec[i] <- vals[2]
  if (i %% 500 == 0) setTxtProgressBar(pb, i)
}
setTxtProgressBar(pb, n); close(pb)

# ==== 결과를 10 m 래스터로 복원 ====
r_ed  <- mask10; terra::values(r_ed)  [cells] <- ed_vec
r_lpi <- mask10; terra::values(r_lpi) [cells] <- lpi_vec
names(r_lpi) <- "LPI100_500m"; names(r_ed) <- "ED100_500m"

# ==== 저장 ====
terra::writeRaster(c(r_lpi, r_ed),
                   file.path(out_dir, "class100_LPI_ED_500m_10m.tif"),
                   overwrite = TRUE, gdal = c("COMPRESS=LZW"))

# (선택) CSV 포인트도 필요하면:
pts <- terra::as.points(mask10, values = FALSE)
pts <- pts[cells]                                     # 계산한 중심만
pts_sf <- sf::st_as_sf(pts)
pts_sf$LPI100_500m <- lpi_vec; pts_sf$ED100_500m <- ed_vec
sf::st_write(pts_sf, file.path(out_dir,"class100_LPI_ED_500m_10m_points.gpkg"),
             delete_dsn = TRUE, quiet = TRUE)

