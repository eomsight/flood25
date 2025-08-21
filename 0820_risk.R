############################
# 0. 공통 세팅
############################

library(terra); library(sf); sf_use_s2(FALSE)

root <- "/Users/hyeseon/Downloads"
crs5179 <- "EPSG:5179"

g30  <- st_read(file.path(root, "grid_30.gpkg"), quiet=TRUE) |> st_transform(5179)
aoi  <- st_as_sf(st_union(st_geometry(g30)))
tmpl <- rast(ext(vect(aoi)), resolution = 30, crs = crs5179)
mask_aoi <- function(r) r |> crop(vect(aoi)) |> mask(vect(aoi))

# 유틸
minmax <- function(r){
  g <- global(r, c("min","max"), na.rm=TRUE)
  out <- (r - g[1,1])/(g[1,2]-g[1,1]); clamp(out, 0, 1, values=TRUE)
}
guess_num_field <- function(sfobj){
  cand <- names(sfobj)[sapply(sfobj, is.numeric)]
  cand <- setdiff(cand, c("fid","id","ID"))
  if (!length(cand)) stop("숫자 필드를 못 찾음"); cand[1]
}
rasterize_to_tmpl <- function(v, field){
  rasterize(vect(v), tmpl, field=field) |> mask_aoi()
}





############################
# 1. Total Vulnerability (기하평균)
############################

# 각 TIF의 마지막 밴드가 V_* 로 저장되어 있다는 전제:
r_soc <- rast(file.path(root, "social_30m.tif"))[["V_social"]]
r_phy <- rast(file.path(root, "physical_30m.tif"))[["V_physical"]]
r_lnd <- rast(file.path(root, "landscape_30m.tif"))[["V_landscape"]]

# 혹시 모를 격자 어긋남 정렬
to30 <- function(r) resample(project(r, crs5179, method="bilinear"), tmpl, method="bilinear") |> mask_aoi()
r_soc <- to30(r_soc); r_phy <- to30(r_phy); r_lnd <- to30(r_lnd)

# 도메인 점수는 0–1 범위가 맞는지 확인 후 보정
r_soc <- minmax(r_soc); r_phy <- minmax(r_phy); r_lnd <- minmax(r_lnd)

V_total <- (r_soc * r_phy * r_lnd)^(1/3)   # 동일가중 기하평균 (SETS 총취약도)
names(V_total) <- "V_total"

writeRaster(c(r_soc, r_phy, r_lnd, V_total),
            file.path(root, "vulnerability_total_30m.tif"),
            overwrite=TRUE, gdal=c("COMPRESS=LZW"))





############################
# 2. Exposure (두 지표 → 30 m 래스터 → 0–1 정규화 → 산술평균)
############################

library(terra); library(sf); sf_use_s2(FALSE)

d <- "/Users/hyeseon/Downloads"
crs5179 <- "EPSG:5179"

g30  <- st_read(file.path(d, "grid_30.gpkg"), quiet=TRUE) |> st_transform(5179)
aoi  <- st_as_sf(st_union(st_geometry(g30)))
tmpl <- rast(ext(vect(aoi)), resolution = 30, crs = crs5179)
mask_aoi <- function(r) r |> crop(vect(aoi)) |> mask(vect(aoi))

minmax <- function(r){
  g <- global(r, c("min","max"), na.rm=TRUE)
  if (is.na(g[1,2] - g[1,1]) || (g[1,2] - g[1,1]) == 0) return(r*0)
  clamp((r - g[1,1])/(g[1,2]-g[1,1]), 0, 1, values=TRUE)
}

# --- 1) 인구밀도: density_num
v_pop <- st_read(file.path(d, "exp_popden.gpkg"), quiet=TRUE) |> st_transform(5179)
if (!"density_num" %in% names(v_pop)) {
  stop(sprintf("exp_popden.gpkg: 'density_num' 필드를 찾을 수 없음. 사용 가능: %s",
               paste(names(v_pop), collapse=", ")))
}
E_pop_raw <- rasterize(vect(v_pop), tmpl, field = "density_num") |> mask_aoi()
E_pop <- minmax(E_pop_raw); names(E_pop) <- "E_popden"

# --- 2) 수계거리: wb_distance  → 근접성(가까울수록 ↑)
v_dist <- st_read(file.path(d, "exp_distance_wb.gpkg"), quiet=TRUE) |> st_transform(5179)
if (!"wb_distance" %in% names(v_dist)) {
  stop(sprintf("exp_distance_wb.gpkg: 'wb_distance' 필드를 찾을 수 없음. 사용 가능: %s",
               paste(names(v_dist), collapse=", ")))
}
E_dist_raw <- rasterize(vect(v_dist), tmpl, field = "wb_distance") |> mask_aoi()
E_prox <- 1 - minmax(E_dist_raw)      # 거리 0→1로 뒤집기(더 가까울수록 값↑)
names(E_prox) <- "E_water_prox"

# --- 3) 노출지수(동일가중 평균) 후 0–1 보정(선택)
E_index <- app(c(E_pop, E_prox), mean, na.rm = TRUE) |> minmax()
names(E_index) <- "E_index"

writeRaster(c(E_pop, E_prox, E_index),
            file.path(d, "exposure_30m.tif"),
            overwrite=TRUE, gdal=c("COMPRESS=LZW"))

# --- 4) 빠른 점검
print(global(c(E_pop, E_prox, E_index), c("min","max","sd"), na.rm=TRUE))








############################ 
############################ 
############################ 여기서부터 다시 ... haz 뭘로 재설정 하느냐에 따라 .. 코드 다시 짜주삼 ㅎㅋ
# 3. Hazard (관악+동작 → 30 m 래스터 → 0/1 정규화)
############################

library(terra); library(sf); sf_use_s2(FALSE)

d <- "/Users/hyeseon/Downloads"
crs5179 <- "EPSG:5179"

# 30 m 템플릿 (grid_30.gpkg 기준)
g30  <- st_read(file.path(d, "grid_30.gpkg"), quiet=TRUE) |> st_transform(5179)
aoi  <- st_as_sf(st_union(st_geometry(g30)))
tmpl <- rast(ext(vect(aoi)), resolution = 30, crs = crs5179)

# 유틸: 0–1 정규화
minmax <- function(r){
  g <- global(r, c("min","max"), na.rm=TRUE)
  clamp((r - g[1,1])/(g[1,2]-g[1,1]), 0, 1, values=TRUE)
}

# 1) 두 구 폴리곤 합치기
hz1 <- st_read(file.path(d, "haz_kwanak.gpkg"),  quiet=TRUE) |> st_make_valid() |> st_transform(5179)
hz2 <- st_read(file.path(d, "haz_dongjak.gpkg"), quiet=TRUE) |> st_make_valid() |> st_transform(5179)
haz <- rbind(hz1, hz2)


# 2) SEG_CODE → 깊이(m) 매핑 (중간값/대표값)
depth_map <- c(N330=0.25, N331=0.75, N332=1.50, N333=3.50, N334=5.00)  # ≥5m는 5로 캡
stopifnot("SEG_CODE" %in% names(haz))
haz$depth_m <- as.numeric(depth_map[haz$SEG_CODE])

# 3) 30 m 래스터화 (폴리곤 밖은 0으로 채움)
H_raw  <- rasterize(vect(haz), tmpl, field = "depth_m") |> crop(vect(aoi)) |> mask(vect(aoi))
H_full <- ifel(is.na(H_raw), 0, H_raw)   # <- 핵심: NA(비침수 지역)를 0으로

# 4) 0–1 정규화
H_norm <- minmax(H_full)

names(H_full) <- "H_depth_m"
names(H_norm) <- "H_norm"

writeRaster(c(H_full, H_norm),
            file.path(d, "hazard_30m.tif"),
            overwrite=TRUE, gdal=c("COMPRESS=LZW"))




############################
# 4. Final Risk = H × E × V_total
############################

# 1) 입력 로드 (AOI 버전 활용)
H_norm   <- rast(file.path(d, "hazard_30m.tif"))[["H_norm"]]
V_total  <- rast(file.path(d, "vulnerability_total_30m.tif"))[["V_total"]]
E_index  <- rast(file.path(d, "exposure_30m_aoi.tif"))[["E_index"]]   # ← 메모리 대신 파일에서

# 2) 정렬(해상도/원점/범위) + AOI 마스크 + 0–1 보정
align <- function(r) resample(project(r, crs5179, method="bilinear"), tmpl, method="bilinear") |> mask_aoi()
H_norm  <- minmax(align(H_norm))
V_total <- minmax(align(V_total))
E_index <- minmax(align(E_index))

# 3) 최종 Risk
Risk <- H_norm * E_index * V_total
Risk <- minmax(Risk); names(Risk) <- "Risk"

# 4) 스택 저장
out_stack <- c(H_norm, E_index, V_total, Risk)
names(out_stack) <- c("H_norm","E_index","V_total","Risk")

writeRaster(out_stack, file.path(d, "risk_30m.tif"),
            overwrite=TRUE, gdal=c("COMPRESS=LZW"))




