### Scalogram: 평균·분산 안정 구간(plateau) 찾기

install.packages("spdep")

library(readr); library(dplyr); library(tidyr); library(stringr); library(ggplot2)
library(sf); library(spdep)

dir <- "/Users/hyeseon/Downloads/outputs"
m   <- read_csv(file.path(dir, "landscape_metrics_by_grid_100m.csv"))
g   <- st_read(file.path(dir, "grid_100m.gpkg"), quiet = TRUE)

# long 포맷으로 변환
mm <- m |>
  pivot_longer(-grid_id, names_to = "name", values_to = "val") |>
  separate(name, into = c("metric","scale"), sep = "_(?=\\d)") |>
  mutate(scale = as.integer(str_remove(scale, "m")))


scalo <- mm %>%
  group_by(metric, scale) %>%
  summarise(mean = mean(val, na.rm=TRUE),
            sd   = sd(val, na.rm=TRUE),
            cv   = sd/mean,
            .groups = "drop") %>%
  arrange(metric, scale)

# 상대변화(%)
scalo_delta <- scalo %>%
  group_by(metric) %>%
  arrange(scale, .by_group=TRUE) %>%
  mutate(rel_change = abs((mean - lag(mean))/lag(mean))*100) %>%
  ungroup()

# 시각화(논문용 그림)
ggplot(scalo, aes(scale, mean)) +
  geom_line() + geom_point() +
  facet_wrap(~metric, scales="free_y") +
  theme_bw() + labs(x="Observation scale (m)", y="Mean metric")

# 간단 룰: 두 구간 모두 2% 미만 -> plateau, 그중 '가장 작은' 스케일 선택
pick_plateau <- scalo_delta %>%
  group_by(metric) %>%
  summarise(
    d1 = rel_change[scale==500],  # 300→500
    d2 = rel_change[scale==700],  # 500→700
    plateau = (!is.na(d1) & !is.na(d2) & d1 < 2 & d2 < 2),
    pick_by_plateau = ifelse(plateau, 300L, NA_integer_)
  )

pick_plateau


### Redundancy: 스케일 간 상관이 너무 높으면 가장 작은 창
cors <- function(pre){ cor(m[grep(paste0("^",pre,"_"), names(m))], use="pairwise") }
cors("SHDI"); cors("LPI"); cors("CONTAG")

