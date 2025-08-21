
# === 필요한 패키지 ===
library(terra); library(sf); sf_use_s2(FALSE)

in_r <- "/Users/hyeseon/Downloads/land_2land.tif"      # 원본(10 m)
in_b <- "/Users/hyeseon/Downloads/boundary.gpkg"       # 경계
out  <- "/Users/hyeseon/Downloads/land_2land_aoi.tif"  # 출력

# 원본 래스터와 경계 읽기
r   <- rast(in_r)
aoi <- st_read(in_b, quiet = TRUE) |>
  st_make_valid() |>
  st_transform(crs(r))            # 래스터 CRS로 맞춤

# 경계로 crop + mask (해상도/원점 보존)
r_aoi <- r |>
  crop(vect(aoi)) |>
  mask(vect(aoi))

# 저장 (NoData 유지, 압축)
writeRaster(r_aoi, out, overwrite = TRUE, gdal = c("COMPRESS=LZW"))

# (선택) 빠른 확인
print(res(r_aoi)); print(ext(r_aoi)); print(crs(r_aoi))
