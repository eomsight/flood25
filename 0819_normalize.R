
#################### social dimension ####################
# ===== 패키지 =====
library(terra); library(sf); sf_use_s2(FALSE)

# ===== 경로/파일 =====
root <- "/Users/hyeseon/Downloads"
f_grid   <- file.path(root, "grid_100.gpkg")
f_senior <- file.path(root, "soc_senior.gpkg")
f_low    <- file.path(root, "soc_low_income.gpkg")
f_dis    <- file.path(root, "soc_disabled.gpkg")
f_shel   <- file.path(root, "soc_dist_shelter.gpkg")

# ===== 템플릿(100 m) & AOI =====
grid <- st_read(f_grid, quiet = TRUE) |> st_transform(5179)
aoi  <- st_as_sf(st_union(st_geometry(grid)))
tmpl <- rast(ext(vect(aoi)), resolution = 100, crs = "EPSG:5179")
mask_aoi <- function(r) r |> crop(vect(aoi)) |> mask(vect(aoi))

# ===== 유틸 =====
guess_num_field <- function(sfobj){
  cand <- names(sfobj)[sapply(sfobj, is.numeric)]
  cand <- setdiff(cand, c("fid","id","ID"))
  if (!length(cand)) stop("숫자 필드를 못 찾음")
  cand[which.min(sapply(sfobj[cand], function(x) sum(is.na(x))))]
}
ras_from_gpkg <- function(gpkg_path, field = NULL){
  v <- st_read(gpkg_path, quiet = TRUE) |> st_transform(5179)
  fld <- if(is.null(field)) guess_num_field(v) else field
  message("Rasterize: ", basename(gpkg_path), " → field='", fld, "'")
  rasterize(vect(v), tmpl, field = fld) |> mask_aoi()
}
minmax <- function(r){
  mn <- global(r, "min", na.rm=TRUE)[1,1]
  mx <- global(r, "max", na.rm=TRUE)[1,1]
  out <- (r - mn)/(mx - mn)
  clamp(out, 0, 1, values = TRUE)
}

# ===== 1) 각 지표: 0–1 정규화 (모두 + 방향) =====
senior   <- ras_from_gpkg(f_senior) |> minmax();    names(senior)   <- "soc_senior"
low_inc  <- ras_from_gpkg(f_low)    |> minmax();    names(low_inc)  <- "soc_low_income"
disabled <- ras_from_gpkg(f_dis)    |> minmax();    names(disabled) <- "soc_disabled"
dist_sh  <- ras_from_gpkg(f_shel)   |> minmax();    names(dist_sh)  <- "soc_dist_shelter"

# 빠른 검증
print(global(c(senior, low_inc, disabled, dist_sh), c("min","max"), na.rm=TRUE))

# ===== 2) Social vulnerability score (산술평균, 동등가중) =====
V_social <- app(c(senior, low_inc, disabled, dist_sh), mean, na.rm = TRUE)
names(V_social) <- "V_social"

# ===== 3) 저장 (QGIS용) =====
stack_out <- c(senior, low_inc, disabled, dist_sh, V_social)
writeRaster(stack_out, file.path(root, "social_100m.tif"),
            overwrite=TRUE, gdal=c("COMPRESS=LZW"))






#################### physical dimension ####################
# ---- 패키지 ----
library(terra); library(sf); sf_use_s2(FALSE)

# ---- 경로 ----
d <- "/Users/hyeseon/Downloads"  # 너의 폴더
g100 <- st_read(file.path(d, "grid_100.gpkg"), quiet = TRUE)
crs5179 <- "EPSG:5179"

# 100 m 템플릿 만들기 (그리드 extent/CRS 기준)
tmpl100 <- rast(ext(vect(g100)), resolution = 100, crs = crs5179)

# ---- 입력 래스터 로드 ----
r_slope <- rast(file.path(d, "phy_slope.tif"))
r_dem   <- rast(file.path(d, "phy_dem.tif"))
r_ndvi  <- rast(file.path(d, "phy_ndvi.tif"))
r_ndbi  <- rast(file.path(d, "phy_ndbi.tif"))

# 전부 5179로 투영 + 100m로 리샘플 (연속형 → bilinear)
to100m <- function(r) {
  if (crs(r) != crs5179) r <- project(r, crs5179, method = "bilinear")
  resample(r, tmpl100, method = "bilinear")
}
r_slope <- to100m(r_slope)
r_dem   <- to100m(r_dem)
r_ndvi  <- to100m(r_ndvi)
r_ndbi  <- to100m(r_ndbi)

# ---- Min–Max 정규화(+뒤집기 옵션) ----
norm01 <- function(x) {
  g <- global(x, c("min","max"), na.rm = TRUE)
  (x - g[1,1]) / (g[1,2] - g[1,1])
}
flip01 <- function(x01) 1 - x01

z_slope <- norm01(r_slope)           # +
z_dem   <- flip01(norm01(r_dem))     # − (뒤집기)
z_ndvi  <- flip01(norm01(r_ndvi))    # − (뒤집기)
z_ndbi  <- norm01(r_ndbi)            # +

names(z_slope) <- "phy_slope"
names(z_dem)   <- "phy_dem_inv"      # 뒤집힌 표고
names(z_ndvi)  <- "phy_ndvi_inv"     # 뒤집힌 NDVI
names(z_ndbi)  <- "phy_ndbi"

# ---- 도메인 점수(가중평균; 기본=동일가중) ----
V_physical <- mean(c(z_slope, z_dem, z_ndvi, z_ndbi), na.rm = TRUE)
names(V_physical) <- "V_physical"

phy <- c(z_slope, z_dem, z_ndvi, z_ndbi, V_physical)

# 경계 마스크가 필요하면 아래 두 줄 주석 해제
aoi <- st_union(st_geometry(g100)) |> st_as_sf() |> st_make_valid()
phy <- mask(phy, vect(aoi))

# ---- 저장 ----
writeRaster(phy, file.path(d, "physical_100m.tif"),
            overwrite = TRUE, gdal = c("COMPRESS=LZW"))



