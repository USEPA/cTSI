# Data Dictionary 

### Overview

This document provides metadata for data used in: Hagy, JD, B. Kreakie, M. Pelletier, F. Nojavan, J. Kiddon, and A. Oczkowski. 2021. Quantifying coastal ecosystem condition and a trophic state index with a Bayesian analytical framework (in review)

The raw data are in the following comma-delimited text files:

[CoastalWQ_20200831.csv](#coastalwq_20200831csv)<br>
[Secchi2015.csv](#Secchi2015.csv)<br>
[bh_secchi.csv](#bh_secchi.csv)<br>
[bh_nutrients.csv](#bh_nutrients.csv)


### CoastalWQ_20200831.csv

This file includes data from EPA's National Coastal Condition Assessments in 2021 and 2015. These data are currently available online at https://www.epa.gov/national-aquatic-resource-surveys/data-national-aquatic-resource-surveys

|	VARIABLE NAME	|	DESCRIPTION 	|	UNITS	|
|:---		|:---		|:---:		|
|	Site_Visnum_Layer 	|	Concatenated unique sample identifier	|		|
|	UID	|	Identifier	|		|
|	SITE_ID           	|	Site name or identifier	|		|
|	VISNUM               	|	Visit number (1 or 2)	|		|
|	Col_Date         	|	Sample collection date	|		|
|	Col_loc          	|	Sample collection location (i.e., surface)	|		|
|	Latitude          	|	Latitude of site	|		|
|	Longitude         	|	Longitude of site	|		|
|	REGION            	|	Region containing site 	|		|
|	SUBREGIONS        	|	Subregion containing site 	|		|
|	STATE             	|	US state containing site	|		|
|	SAMPYEAR         	|	Year sample was collected	|		|
|	CHLA..ug.L.           	|	Chlorophyll-a 	|	micrograms per liter	|
|	TN..mgN.L.        	|	Total nitrogen 	|	milligrams N per liter	|
|	DIN..mgN.L.      	|	Dissolved inorganic nitrogen 	|	milligrams N per liter	|
|	TP..mgP.L.       	|	Total phosphorus 	|	milligrams P per liter	|
|	DIP..mgP.L.       	|	Dissolved inorganic phosphorus 	|	milligrams P per liter	|
|	DIN.DIP..Molar.  	|	Molar ratio of dissolved inorganic N to P	|		|
|	TN.TP..Molar.     	|	Molar ratio of total N to P	|		|
|	SECCHI_MEAN..m.       	|	Average of secchi depth measurements during the visit 	|	meters	|
|	TSS..mg.L.        	|	Total suspended solids 	|	milligrams per liter	|
|	levelIII              	|	Level III ecoregion containing site	|		|
|	AggEco                	|	Aggregate ecoregion containing site	|		|
|	EMAP             	|	Environmental Monitoring and Assessment Program region containing site	|		|
|	BioGeo	|	Biogeographic province containing site	|		|


### Secchi2015.csv

Includes secchi depth measurements from the 2015 NCCA assessment that were not included in the initial file provided to the investigators.

|	VARIABLE NAME	|	DESCRIPTION 	|	UNITS	|
|:---		|:---		|:---:		|
|SITE_ID	|		Site name or identifier | |
|VISIT_NO	|	Visit number (1 or 2) | |
|MEAN_SECCHI_DEPTH |   Average of secchi depth measurements | meters |
|Clear to bottom	|	Flag indicating "Y" is secchi disk was visible on the bottom| |


### bh_secchi.csv

Water clarity data from Boston Harbor downloaded from MWRA at https://www.mwra.com/harbor/html/wq_data.htm

|	VARIABLE NAME	|	DESCRIPTION 	|	UNITS	|
|:---		|:---		|:---:		|
|	Project.ID	|	MWRA project identifier	|		|
|	Region                                  	|	MWRA region containing station	|		|
|	Subregion                               	|	MWRA sub-region containing station	|		|
|	DEP.Segment                              	|	Massachusetts DEP segment identifier	|		|
|	Station.ID                               	|	Station number	|		|
|	Surface.or.Bottom                       	|	depth layer of samples (S= Surface, B= Bottom)	|		|
|	Depth.to.bottom..m.                      	|	Water depth	|	meters	|
|	Sample.depth.for.TSS..m.                 	|	Depth where total suspended solids sample was collected	|		|
|	Date.time..EASTERN.STANDARD.TIME.       	|	Sample time in Eastern Standard Time	|		|
|	Total.Suspended.Solids..mg.L.            	|	Total suspended solids concentration	|	milligrams per liter	|
|	Total.Suspended.Solids............or.... 	|	Flag for total suspended solids out of range (blank, > or <)	|		|
|	Secchi.disk.depth..m.                   	|	Secchi disk depth	|	meters	|
|	Attenuation.Coefficient..m.1.            	|	Light attentuation coefficients	|	per meter	|

### bh_nutrients.csv

Nutrients data from Boston Harbor downloaded from MWRA at https://www.mwra.com/harbor/html/wq_data.htm


|	Project.ID							|	MWRA project identifier	|		|
|:---		|:---		|:---:		|
|	Region	|			MWRA region containing station	|		|
|	Subregion	|			MWRA sub-region containing station	|		|
|	DEP.Segment	|			Massachusetts DEP segment identifier	|		|
|	Station.ID	|			Station number	|		|
|	Date.time..EASTERN.STANDARD.TIME.			|		Time of sample collected in Eastern Standard Time	|		|
|	Surface.or.Bottom	|			depth layer of samples (S= Surface, B= Bottom)	|		|
|	Sample.depth..m.                  				|	Depth where sample was collected (meters)	|	meters	|
|	Ammonium..uM.                    				|	Ammonium concentration	|	micromoles per liter	|
|	Ammonium.....or....               				|	Ammonium below detection limit (<)	|		|
|	Nitrate.nitrite..uM.              		|			Nitrate plus nitrite	|	micromoles per liter	|
|	Nitrate.nitrite......or....      		|			Nitrate plus nitrite below detection limit (<)	|		|
|	Total.dissolved.N..uM.            		|			Total dissolved nitrogen 	|	micromoles N per liter	|
|	Total.dissolved.N......or....	|				Total dissolved nitrogen below detection limit (<)	|		|
|	Particulate.nitrogen..uM.       		|			Particulate nitrogen	|	micromoles N per liter	|
|	Particulate.nitrogen......or....	|				Particulate nitrogen below detection limit (<)	|		|
|	Total.Kjeldahl.Nitrogen..uM.     		|			Total Kjeldahl nitrogen 	|	micromoles N per liter	|
|	TKN............or....	|				Total Kjeldahl nitrogen below detection limit (<)	|		|
|	Phosphate..uM.	|				Phosphate	|	micromoles P per liter	|
|	Phosphate.....or....	|				Phosphate below detection limit (<)	|		|
|	Total.dissolved.P..uM.					|		Total dissolved phosphorus	|	micromoles P per liter	|
|	Total.dissolved.P.....or....      		|			Total dissolved phosphorus below detection limit (<)	|		|
|	Particulate.P..uM.	|				Particulate phosphorus	|	micromoles P per liter	|
|	Particulate.P.....or....         		|			Particulate phosphorus below detection limit (<)	|		|
|	Total.phosphorus..uM.	|				Total phosphorus	|	micromoles P per liter	|
|	Total.phosphorus.........or....	|				Total phosphorus below detection limit (<)	|		|
|	Particulate.carbon..uM.          		|			Particulate carbon	|	micromoles C per liter	|
|	Particulate.carbon.........or....	|				Particulate carbon below detection limit (<)	|		|
|	Chlorophyll.a..ug.L.	|				Chlorophyll-a	|	micromoles per liter	|
|	Chlorophyll.....or....	|				Chlorophyll-a below detection limit (<)	|		|
|	Phaeophytin..ug.L.	|				Phaeophytin	|	micromoles per liter	|
|	Phaeophytin.....or....					|		Phaeophytin below detection limit (<)	|		|

