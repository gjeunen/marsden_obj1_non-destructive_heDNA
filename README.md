# Marsden objective 1: non-destructive heDNA

## 1. Introduction

This GitHub repository serves as a comprehensive guide for the bioinformatic and statistical analysis of the ethanol comparison project, associated with the Marsden Fast-Start funding MFP-UOO002116.

## 2. Experimental design

The primary objective of this experiment is to assess the possibility of non-destructive heDNA recovery from museum-stored marine sponge specimens through DNA extraction from the ethanol in which specimens are stored. To verify non-destructive sampling as a suitable methodological choice for specimens across the phylum Porifera, we selected one specimen from the three major classes due to differing internal skeletal structures, including Hexactinellida (hexagonal silica spicules), Demospongiae (spongin), and Calcarea (calcium carbonate spicules). A detailed list with metadata for each specimen can be found below.

| Cat. # | Class | Species | Date | Latitude | Longitude | Depth |
| :---- | :---- | :---- | :---- | :---- | :---- | :---- |
| 35574 | Calcarea | Petrobiona sp. | 09/02/2008 | -73.1245 | 174.3205 | 321 |
| 37535 | Demospongia | Homaxinella sp. | 24/02/2008 | -72.0093 | 173.2238 | 850 |
| 37124 | Hexactinellida | Rossella fibulata | 21/02/2008 | -72.5903 | 175.3423 | 479 |

For each sepcimen, DNA was extracted from 10 tissue biopsies without removal of excess ethanol and 10 tissue biopsies whereby excess ethanol was removed by absorbance of lint-free kimwipes. To test heDNA recovery from ethanol, we extracted DNA through (1) 1 ml evaporation, (2) 1 ml centrifugation, (3) 1 ml precipitation, (4) 1 ml filtration, and (5) 10 ml filtration. For each non-destructive method, 10 replicates were processed to enable statistical comparison between treatments.

After DNA extraction, the total DNA concentration was measured via Qubit, while DNA purity was assessed by investigating the 260/280 and 260/230 absorbence ratios measured via Denovix. To amplify heDNA signals associated with fish species (Fish16SF: 5’-GACCCTATGGAGCTTTAGAC-3’; Fish16S2R: 5’-CGCTGTTATCCCTADRGTAACT-3’), we employed a metabarcoding approach. This primer set was chosen over the MiFish-U primer set, as it does not co-amplify human DNA, a potential contaminant signal for museum-stored specimens that have been handled repeatedly for morphometric analysis.

**While we're unsure at this stage, genome sequencing was also conducted on the ONT platform to assess the possibility of retrieving (mito) genome sequences from non-destructive methods. However, we'll have to wait until we receive the results later on if this will be included in the manuscript.**

## 3. Figure 1: Map of Antarctica

The first figure of the manuscript is a map of Antarctica displaying the locations of the three specimens which were analysed in this experiment. To generate the map, we can run the R script below. To successfully execute the script, several files will need to be downloaded, including:

1. A high resolution Antarctic outline including ice sheets:
   1. source: <https://data.bas.ac.uk/full-record.php?id=GB/NERC/BAS/PDC/01391>
   2. filename: add_coastline_high_res_polygon_v7.2.shp
2. A csv file containing the specimen metadata
   1. filename: specimen_metadata.csv
3. A high resolution bathymetry data file as baselayer:
   1. source: <https://www.gebco.net>
   2. filename: IBCSO_v2_ice-surface.nc
4. A world map with borders for all countries:
   1. source: <https://www.star.nesdis.noaa.gov/data/smcd1/vhp/GIS/TM_World_Borders/>
   2. filename: TM_WORLD_BORDERS-0.3.shp

All files needed to execute the R code can be found in the subdirectory `metadata/figures/figure_1/` within this GitHub repo. Once the figure is generated in R and exported to a pdf file, import figure into Adobe Illustrator to clean up figure for publication.

```{code-block} R
#######################
## MAP OF ANTARCTICA ##
#######################
# set working directory and load libraries
setwd("")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(terra, tidyterra, ggplot2, ggnewscale, patchwork)
require(sf)
library(recolorize)
library(egg)
library(png)
library(ggpubr)
library(viridis)
library(RColorBrewer) 
library(dplyr)
library(hrbrthemes)
library(readxl)
library(ggrepel)
library(pracma)
library(scales)

# read data
ATAshp <- simplifyGeom(vect("add_coastline_high_res_polygon_v7.2/add_coastline_high_res_polygon_v7.2.shp"), tolerance=500)
ATAshp$col <- sapply(ATAshp$surface, function(x){ifelse(x == "land", "grey85", "grey98")})
sponge_dat <- read.csv("specimen_metadata2.csv")
sponge_pts <- project(vect(sponge_dat, geom = c("Longitude1", "Latitude1"), keepgeom = TRUE, crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"), ATAshp)
grat <- project(vect(sf::st_graticule(lon=seq(-175, 180, 5), lat = seq(-85, -60, 5), ndiscr = 5000)), ATAshp)
ibsco <- rast("add_coastline_high_res_polygon_v7.2/IBCSO_v2_ice-surface.nc")
sample_colors <- c("Calcarea" = "lightgoldenrod", "Demospongiae" = "steelblue", "Hexactinellida" = "firebrick")
spp_names <- c(expression(italic("Calcarea")), expression(italic("Demospongiae")), expression(italic("Hexactinellida")))
sponge_pts$Catalog.Number <- as.factor(sponge_pts$Catalog.Number)
xmn <- -1000000
xmx <- 1500000
ymn <- -2500000
ymx <- -1000000
WORLDshp <- simplifyGeom(vect("TM_WORLD_BORDERS-0.3/TM_WORLD_BORDERS-0.3.shp"))
WORLDshpNoATA <- WORLDshp[WORLDshp$ISO3 != "ATA"]

# generate circular masking layer for overview plot
prj <- "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs +ellps=WGS84 +towgs84=0,0,0"
little <- st_point(c(0,0)) %>% st_sfc(crs = prj) %>% st_buffer(dist = 4530)
xlim <- sf::st_bbox(little)[c("xmin", "xmax")]*2.5
ylim <- sf::st_bbox(little)[c("ymin", "ymax")]*2.5
encl_rect <- list(cbind(c(xlim[1], xlim[2], xlim[2], xlim[1], xlim[1]), 
                        c(ylim[1], ylim[1], ylim[2], ylim[2], ylim[1]))) %>%
  sf::st_polygon() %>%
  sf::st_sfc(crs = prj)
circleMask <- sf::st_difference(encl_rect, little)

# generate plots
RSR_plot <- ggplot() + 
  geom_spatraster(data = ibsco, show.legend = FALSE) + 
  scale_fill_gradient(low = "#528B8B", high = "#CBDCDC") + 
  geom_spatvector(data = grat, col = "grey50", alpha = 0.2) + 
  geom_spatvector(data = ATAshp, alpha = 1, fill = ATAshp$col, col = 'grey20') + 
  new_scale_fill() + 
  geom_spatvector(data = sponge_pts, aes(shape = Class, fill = Class), alpha = 1, size = 5, stroke = 0.5) + 
  scale_fill_manual(values = sample_colors, labels = c("Calcarea", "Demospongiae", "Hexactinellida") , name = "Sponge Class") +
  guides(fill = guide_legend(override.aes = list(pch=21))) + 
  scale_shape_manual(values = c(21, 22, 23)) + 
  scale_x_continuous(breaks = seq(-180, 180, by = 10)) + 
  coord_sf(crs = crs(ATAshp), expand = FALSE, xlim = c(xmn, xmx), ylim = c(ymn, ymx)) + 
  annotate(geom = "rect", xmin = xmn, xmax = xmx, ymin = ymn, ymax = ymx, fill = NA, col = "grey10") +  
  xlab(NULL) + ylab(NULL) #  

inset_map <- ggplot() + 
  geom_spatvector(data = ATAshp, alpha = 1, fill = ATAshp$col, col = 'grey20') + 
  geom_rect(aes(xmin = xmn, xmax = xmx, ymin = ymn, ymax = ymx), fill = "dark red", col = "black", alpha = 0.2) + 
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), plot.background = element_blank(), panel.border = element_blank(), axis.line = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_blank()) 

#output to pdf
pdf("RSR_spongemap_v2.pdf", width = 8, height = 5.5) 
RSR_plot + annotation_custom(ggplotGrob(inset_map), xmin = xmx - (xmx-xmn)/3, xmax = xmx, ymin = ymx - (ymx-ymn)/3, ymax = ymx + 280000) 
dev.off() 
```
