# Definitions
# soil_var <- "clay"
soil_var <- "phh2o"
if (soil_var == "clay") {
  soil_var_color <- RColorBrewer::brewer.pal(9, "Oranges")
  soil_var_unit <- "g/kg"
  soil_var_name <- "Argila"
} else if (soil_var == "phh2o") {
  soil_var_color <- RColorBrewer::brewer.pal(9, "Blues")
  soil_var_unit <- ". * 10"
  soil_var_name <- "pH"
}
which_biome <- "Caatinga"
# which_biome <- "Cerrado"
# which_biome <- "Pantanal"

# Download vector of Brazilian Biomes
if (!require(geobr)) {
  install.packages("geobr")
}
biome <- geobr::read_biomes()
biome <- biome[biome[["name_biome"]] == which_biome, "name_biome"]
sf::st_write(biome, "data/cutline.geojson", delete_dsn = TRUE)
plot(biome)

# FEBR #############################################################################################
if (!require(data.table)) {
  install.packages("data.table")
}
febr <- data.table::fread(
  "/home/alessandrorosa/ownCloud/febr-repo/publico/febr-superconjunto.txt", dec = ",")
febr <- febr[profund_sup < 20, ]
febr <- febr[!is.na(coord_x), ]
febr <- febr[!is.na(coord_y), ]
febr[, profund_inf := ifelse(profund_inf > 20, 20, profund_inf)]
febr[, thickness := profund_inf - profund_sup]
febr[, id := paste0(dataset_id, "-", observacao_id)]
febr_cols <- c("id",
  # "profund_sup", "profund_inf",
  "thickness", "coord_x", "coord_y",
  # "ctc", "carbono", "silte", "areia",
  "argila", "ph", "dsi")
febr <- febr[, ..febr_cols]
febr[, clay := argila * thickness / 20]
febr[, phh2o := ph * thickness / 20]
febr[, phh2o := ph * 10]
febr[, bdod := dsi * thickness / 20]
febr <- febr[, c("id", "coord_x", "coord_y", "clay", "phh2o", "bdod")]
febr_mean <- febr[, lapply(.SD, sum, na.rm = TRUE), by = c("id", "coord_x", "coord_y")]
febr_mean[, clay := ifelse(clay < 10, NA_real_, clay)]
febr_mean[, clay := ifelse(clay > 1000, NA_real_, clay)]
febr_mean[, bdod := ifelse(bdod < 0.5, NA_real_, bdod)]
febr_mean[, bdod := ifelse(bdod > 2, NA_real_, bdod)]
febr_mean[, phh2o := ifelse(phh2o < 20, NA_real_, phh2o)]
febr_mean[, phh2o := ifelse(phh2o > 80, NA_real_, phh2o)]
# count FEBR data points
length2 <- function(x) {length(na.exclude(x))}
febr_mean[, lapply(.SD, length2)]
# plot FEBR points
febr_mean <- sf::st_as_sf(febr_mean, coords = c("coord_x", "coord_y"), crs = 4674)
# dev.off()
png("res/febr-points.png", width = 480 * 2, height = 480 * 2, res = 72 * 2)
plot(febr_mean, graticule = TRUE, axes = TRUE)
dev.off()

# SOILGRIDS ########################################################################################
# https://git.wur.nl/isric/soilgrids/soilgrids.notebooks/-/blob/master/markdown/webdav_from_R.md
# https://files.isric.org/soilgrids/latest/data/
homolosine <- "+proj=igh +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs"
biomes_homolosine <- sf::st_transform(biome, crs = homolosine)
soilgrids_depths <- c("0-5", "5-15", "15-30")
soil_files <- paste0(soil_var, "/", soil_var, "_", soilgrids_depths, "cm_mean.vrt")
sg_url <- paste0(
  "/vsicurl?max_retry=3&retry_delay=1&list_dir=no&url=",
  "https://files.isric.org/soilgrids/latest/data/", soil_files)
for (i in sg_url) {
  gdalUtils::gdal_translate(
    src_dataset = i,
    dst_dataset = paste0("tmp/", basename(i)),
    tr = c(250, 250), of = "VRT",
    projwin = sf::st_bbox(biomes_homolosine)[c(1, 4, 3, 2)],
    projwin_srs = homolosine, verbose = TRUE)
  gdalUtils::gdalwarp(
    srcfile = paste0("tmp/", basename(i)),
    dstfile = paste0("tmp/warped-", basename(i)),
    of = "VRT",
    s_srs = homolosine,
    t_srs = "EPSG:4674")
  gdalUtils::gdal_translate(
    src_dataset = paste0("tmp/warped-", basename(i)),
    of = "GTiff",
    dst_dataset = paste0("data/", strsplit(basename(i), ".", fixed = TRUE)[[1]][1], ".tif"))
}
filename <- paste0("data/", basename(sg_url))
filename <- gsub(".vrt", ".tif", filename, fixed = TRUE)
filename <- paste0("-", LETTERS[1:3], " ", filename, collapse = " ")
outfile <- paste0("data/", tolower(which_biome), "_", soil_var, "_0-20cm_mean.tif")
cmd <- paste0("gdal_calc.py --calc A * 0.25 + B * 0.5 + C * 0.25 ", filename, " --outfile ", outfile)
system(cmd)
outfile_cropped <- gsub("data/", "data/cropped_", outfile)
cmd <- paste0("gdalwarp -cutline data/cutline.geojson -crop_to_cutline -co COMPRESS=DEFLATE ",
  outfile, " ", outfile_cropped)
system(cmd)
# Plot SoilGrids
soilgrids <- stars::read_stars(outfile_cropped)
febr_mean_in <- sf::st_intersects(biome, febr_mean, sparse = FALSE)
febr_mean <- febr_mean[febr_mean_in, ]
breaks <- range(soilgrids[[1]], na.rm = TRUE)
breaks <- seq(breaks[1], breaks[2], length.out = 10)
png_file <- gsub(".tif", ".png", basename(outfile_cropped), fixed = TRUE)
dev.off()
png(paste0("res/", png_file), width = 480 * 2, height = 480 * 2, res = 72 * 2)
plot(soilgrids, axes = TRUE, col = soil_var_color, breaks = breaks,
  reset = FALSE, main = "")
plot(febr_mean[soil_var], add = TRUE, col = "black")
dev.off()
# Validate SoilGrids
febr_mean$soilgrids <- stars::st_extract(soilgrids, febr_mean)[[1]]
error <- na.exclude(febr_mean[["soilgrids"]] - febr_mean[[soil_var]])
residual <- na.exclude(
  mean(febr_mean[[soil_var]], na.rm = TRUE) - febr_mean[[soil_var]])
me <- round(mean(error))
rmse <- round(sqrt(mean(sum(error * error))))
nse <- 1 - round(sum(error * error) / sum(residual * residual), 2)
xlim <- ylim <- range(c(febr_mean[[soil_var]], febr_mean[["soilgrids"]]), na.rm = TRUE)
dev.off()
png_file <- paste0("res/", tolower(soil_var_name), "_scatterplot.png")
png(png_file, width = 480 * 2, height = 480 * 2, res = 72 * 2)
plot(x = febr_mean[["soilgrids"]],
  y = febr_mean[[soil_var]],
  ylim = ylim, xlim = xlim, panel.first = grid(),
  ylab = paste0(soil_var_name, " (", soil_var_unit, ")"),
  xlab = paste0("SoilGrids250m (", soil_var_unit, ")"))
if (soil_var == "phh2o") {
  mtext(text = paste0(soil_var_name, " na camada 0-20 cm"),
    adj = 0, at = 35, cex = 1.5, line = 2, font = 2)
  mtext(text = "Ajuste entre SoilGrids250m v2.0 e dados observados obtidos do FEBR",
    adj = 0, at = 35, cex = 1, line = 0.5, font = 1)
  text(x = xlim[2] * 0.85, y = ylim[2] * 0.6, pos = 4,
    label = paste0("ME = ", me, " ", soil_var_unit))
  text(x = xlim[2] * 0.85, y = ylim[2] * 0.575, pos = 4,
    label = paste0("RMSE = ", rmse, " ", soil_var_unit))
  text(x = xlim[2] * 0.85, y = ylim[2] * 0.55, pos = 4,
    label = paste0("NSE = ", nse))
} else {
  mtext(text = paste0("Conteúdo de ", tolower(soil_var_name), " na camada 0-20 cm"),
    adj = 0, at = 24, cex = 1.5, line = 2, font = 2)
  mtext(text = "Ajuste entre SoilGrids250m v2.0 e dados observados obtidos do FEBR",
    adj = 0, at = 24, cex = 1, line = 0.5, font = 1)
  text(x = xlim[2] * 0.7, y = ylim[2] * 0.2, pos = 4,
    label = paste0("ME = ", me, " ", soil_var_unit))
  text(x = xlim[2] * 0.7, y = ylim[2] * 0.15, pos = 4,
    label = paste0("RMSE = ", rmse, " ", soil_var_unit))
  text(x = xlim[2] * 0.7, y = ylim[2] * 0.1, pos = 4,
    label = paste0("NSE = ", nse))
}
abline(a = 0, b = 1, col = "blue", lty = "dashed")
fit <- lm(as.formula(paste0(soil_var, "~ soilgrids")), febr_mean)
abline(fit, col = "red", lty = "dashed")
dev.off()
# # Geostatistical interpolation
# febr_mean$residual <- febr_mean[["soilgrids"]] - febr_mean[[soil_var]]
# is_na_residual <- is.na(febr_mean$residual)
# febr_mean_vario <- febr_mean[!is_na_residual, ]
# febr_vario <- gstat::variogram(residual ~ 1, sf::as_Spatial(febr_mean_vario))
# febr_vario_model <- gstat::vgm(model = "Exp", psill = 2000, nugget = 11000, range = 100)
# plot(febr_vario, model = febr_vario_model)
# tmp <- gstat::krige(residual ~ 1,
#   sf::as_Spatial(febr_mean_vario),
#   as(as(soilgrids, "Raster"), "SpatialPointsDataFrame"),
#   model = febr_vario_model)

