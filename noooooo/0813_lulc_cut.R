install.packages("lwgeom")

library(sf)
library(lwgeom)

# 경로 지정
lulc_path     <- "/Users/hyeseon/Downloads/lulc_1m_whole.shp"
boundary_path <- "/Users/hyeseon/Downloads/boundary.gpkg"
out_path      <- "/Users/hyeseon/Downloads/LULC_boundary.gpkg"

# 데이터 읽기
lulc <- st_read(lulc_path, quiet = TRUE)
boundary <- st_read(boundary_path, quiet = TRUE)

# 좌표계 맞추기
if (st_crs(lulc) != st_crs(boundary)) {
  boundary <- st_transform(boundary, st_crs(lulc))
}

# 지오메트리 오류 수정
lulc <- st_make_valid(lulc)
boundary <- st_make_valid(boundary)

# 경계 dissolve
boundary_union <- st_union(boundary)

# 클립
lulc_clipped <- st_intersection(lulc, boundary_union)

# 저장 (gpkg)
if (file.exists(out_path)) file.remove(out_path)
st_write(lulc_clipped, out_path, layer = "LULC_boundary", quiet = TRUE)

cat("저장 완료:", out_path, "\n")
