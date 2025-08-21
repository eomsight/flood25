
#################### social dimension (30 m) ####################
# ===== 패키지 =====
library(terra); library(sf); sf_use_s2(FALSE)

# ===== 경로/파일 =====
root <- "/Users/hyeseon/Downloads"
f_grid   <- file.path(root, "grid_30.gpkg")              # <<< 30 m 그리드
f_senior <- file.path(root, "soc_senior.gpkg")
f_low    <- file.path(root, "soc_low_income.gpkg")
f_dis    <- file.path(root, "soc_disabled.gpkg")
f_shel   <- file.path(root, "soc_dist_shelter.gpkg")     # 거리: 멀수록 취약 ↑ (양의 관계)

# ===== 템플릿(30 m) & AOI =====
grid <- st_read(f_grid, quiet = TRUE) |> st_transform(5179)
aoi  <- st_as_sf(st_union(st_geometry(grid)))
tmpl <- rast(ext(vect(aoi)), resolution = 30, crs = "EPSG:5179")   # <<< 30 m 템플릿
mask_aoi <- function(r) r |> crop(vect(aoi)) |> mask(vect(aoi))

# ===== 유틸 =====
ras_from_gpkg <- function(gpkg_path, field){
  v <- st_read(gpkg_path, quiet = TRUE) |> st_transform(5179)
  if (!field %in% names(v)) {
    stop(sprintf("%s: 지정한 필드 '%s'가 없음. 사용 가능: %s",
                 basename(gpkg_path), field, paste(names(v), collapse=", ")))
  }
  message("Rasterize: ", basename(gpkg_path), " → field='", field, "'")
  rasterize(vect(v), tmpl, field = field) |> mask_aoi()
}
minmax <- function(r){
  mn <- global(r, "min", na.rm=TRUE)[1,1]
  mx <- global(r, "max", na.rm=TRUE)[1,1]
  if (is.na(mx - mn) || (mx - mn) == 0) return(r*0)
  out <- (r - mn)/(mx - mn)
  clamp(out, 0, 1, values = TRUE)
}

# ===== 1) 각 지표: 0–1 정규화 (모두 + 방향) =====
senior   <- ras_from_gpkg(f_senior, "X65_percent_num")     |> minmax(); names(senior)   <- "soc_senior"
low_inc  <- ras_from_gpkg(f_low,    "low_percent_num")     |> minmax(); names(low_inc)  <- "soc_low_income"
disabled <- ras_from_gpkg(f_dis,    "disabled_percent_num")|> minmax(); names(disabled) <- "soc_disabled"
dist_sh  <- ras_from_gpkg(f_shel,   "shelter_Distance")    |> minmax(); names(dist_sh)  <- "soc_dist_shelter"

# 빠른 검증
print(global(c(senior, low_inc, disabled, dist_sh), c("min","max"), na.rm=TRUE))

# ===== 2) Social vulnerability score (산술평균, 동등가중) =====
V_social <- app(c(senior, low_inc, disabled, dist_sh), mean, na.rm = TRUE)
names(V_social) <- "V_social"

# ===== 3) 저장 (QGIS용 GeoTIFF) =====
stack_out <- c(senior, low_inc, disabled, dist_sh, V_social)
writeRaster(stack_out, file.path(root, "social_30m.tif"),
            overwrite=TRUE, gdal=c("COMPRESS=LZW"))






#################### physical dimension (30 m) ####################
# ---- 패키지 ----
library(terra); library(sf); sf_use_s2(FALSE)

# ---- 경로 ----
d <- "/Users/hyeseon/Downloads"
g30 <- st_read(file.path(d, "grid_30.gpkg"), quiet = TRUE)
crs5179 <- "EPSG:5179"

# 30 m 템플릿 (그리드 extent/CRS 기준)
tmpl30 <- rast(ext(vect(g30)), resolution = 30, crs = crs5179)

# ---- 입력 래스터 로드 ----
r_slope <- rast(file.path(d, "phy_slope.tif"))
r_dem   <- rast(file.path(d, "phy_dem.tif"))
r_ndvi  <- rast(file.path(d, "phy_ndvi.tif"))
r_ndbi  <- rast(file.path(d, "phy_ndbi.tif"))

# 전부 5179로 투영 + 30m로 리샘플 (연속형 → bilinear)
to30m <- function(r) {
  if (crs(r) != crs5179) r <- project(r, crs5179, method = "bilinear")
  resample(r, tmpl30, method = "bilinear")
}
r_slope <- to30m(r_slope)
r_dem   <- to30m(r_dem)
r_ndvi  <- to30m(r_ndvi)
r_ndbi  <- to30m(r_ndbi)

# ---- Min–Max 정규화(+뒤집기 옵션) ----
norm01 <- function(x) {
  g <- global(x, c("min","max"), na.rm = TRUE)
  out <- (x - g[1,1]) / (g[1,2] - g[1,1])
  clamp(out, 0, 1, values = TRUE)
}
flip01 <- function(x01) 1 - x01

# 부호관계: slope(+) / dem(-) / ndvi(-) / ndbi(+)
z_slope <- norm01(r_slope)            # +
z_dem   <- flip01(norm01(r_dem))      # −
z_ndvi  <- flip01(norm01(r_ndvi))     # −
z_ndbi  <- norm01(r_ndbi)             # +

names(z_slope) <- "phy_slope"
names(z_dem)   <- "phy_dem_inv"       # 뒤집힌 DEM
names(z_ndvi)  <- "phy_ndvi_inv"      # 뒤집힌 NDVI
names(z_ndbi)  <- "phy_ndbi"

# ---- 도메인 점수(동일가중 산술평균) ----
V_physical <- mean(c(z_slope, z_dem, z_ndvi, z_ndbi), na.rm = TRUE)
names(V_physical) <- "V_physical"

phy <- c(z_slope, z_dem, z_ndvi, z_ndbi, V_physical)

# 경계 마스크(그리드 union)
aoi <- st_union(st_geometry(g30)) |> st_as_sf() |> st_make_valid()
phy <- mask(phy, vect(aoi))

# ---- 저장 ----
writeRaster(phy, file.path(d, "physical_30m.tif"),
            overwrite = TRUE, gdal = c("COMPRESS=LZW"))






#################### landscape dimension (30 m) ####################
library(terra); library(sf); sf_use_s2(FALSE)

# ---- 경로/그리드/템플릿 ----
d <- "/Users/hyeseon/Downloads"
g30 <- st_read(file.path(d, "grid_30.gpkg"), quiet = TRUE)
crs5179 <- "EPSG:5179"
tmpl30 <- rast(ext(vect(g30)), resolution = 30, crs = crs5179)

# ---- 입력 래스터 로드 (필수 밴드만 선택) ----
r_land  <- rast(file.path(d, "land_2land_aoi.tif"))    # band1=SHDI, band3=CONTAG
r_class <- rast(file.path(d, "land_2class.tif"))   # band1=LPI,  band2=ED

r_shdi   <- r_land[[1]]
r_contag <- r_land[[3]]
r_lpi    <- r_class[[1]]
r_ed     <- r_class[[2]]

# ---- 5179 투영 + 30 m 리샘플 ----
to30m <- function(r){
  if (crs(r) != crs5179) r <- project(r, crs5179, method = "bilinear")
  resample(r, tmpl30, method = "bilinear")
}
r_shdi   <- to30m(r_shdi)
r_contag <- to30m(r_contag)
r_lpi    <- to30m(r_lpi)
r_ed     <- to30m(r_ed)

# ---- 0–1 정규화(+/− 반영) ----
norm01 <- function(x){
  g <- global(x, c("min","max"), na.rm = TRUE)
  out <- (x - g[1,1]) / (g[1,2] - g[1,1])
  clamp(out, 0, 1, values = TRUE)
}
flip01 <- function(x01) 1 - x01


z_shdi   <- flip01(norm01(r_shdi))      # − (SHDI↑ ⇒ 취약↓)
z_ed     <- flip01(norm01(r_ed))        # − (ED↑   ⇒ 취약↓)
z_contag <- norm01(r_contag)            # + (CONTAG↑ ⇒ 취약↑)
z_lpi    <- norm01(r_lpi)               # + (LPI↑    ⇒ 취약↑)

names(z_shdi)   <- "land_shdi_inv"
names(z_ed)     <- "land_ed_inv"
names(z_contag) <- "land_contag"
names(z_lpi)    <- "land_lpi"


# ---- 도메인 점수(동일가중 산술평균) ----
V_landscape <- mean(c(z_shdi, z_ed, z_contag, z_lpi), na.rm = TRUE)
names(V_landscape) <- "V_landscape"

land <- c(z_shdi, z_ed, z_contag, z_lpi, V_landscape)

# 그리드 경계로 마스크
aoi <- st_union(st_geometry(g30)) |> st_as_sf() |> st_make_valid()
land <- mask(land, vect(aoi))

# ---- 저장 ----
writeRaster(land, file.path(d, "landscape_30m.tif"),
            overwrite = TRUE, gdal = c("COMPRESS=LZW"))

# 빠른 검증(선택)
print(global(land, c("min","max"), na.rm = TRUE))
