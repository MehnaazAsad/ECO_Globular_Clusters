ECO_Globular_Clusters
==============================

Repository for ECO project regarding Globular Clusters

Project Organization
------------

    ├── LICENSE
    ├── Makefile           <- Makefile with commands like `make data` or `make train`
    ├── README.md          <- The top-level README for developers using this project.
    ├── data
    │   ├── external       <- Data from third party sources.
    │   ├── interim        <- Intermediate data that has been transformed.
    │   ├── processed      <- The final, canonical data sets for modeling.
    │   └── raw            <- The original, immutable data dump.
    │
    ├── docs               <- A default Sphinx project; see sphinx-doc.org for details
    │
    ├── models             <- Trained and serialized models, model predictions, or model summaries
    │
    ├── notebooks          <- Jupyter notebooks. Naming convention is a number (for ordering),
    │                         the creator's initials, and a short `-` delimited description, e.g.
    │                         `1.0-jqp-initial-data-exploration`.
    │
    ├── references         <- Data dictionaries, manuals, and all other explanatory materials.
    │
    ├── reports            <- Generated analysis as HTML, PDF, LaTeX, etc.
    │   └── figures        <- Generated graphics and figures to be used in reporting
    │
    ├── requirements.txt   <- The requirements file for reproducing the analysis environment, e.g.
    │                         generated with `pip freeze > requirements.txt`
    │
    ├── src                <- Source code for use in this project.
    │   ├── __init__.py    <- Makes src a Python module
    │   │
    │   ├── data           <- Scripts to download or generate data
    │   │   └── make_dataset.py
    │   │
    │   ├── features       <- Scripts to turn raw data into features for modeling
    │   │   └── build_features.py
    │   │
    │   ├── models         <- Scripts to train models and then use trained models to make
    │   │   │                 predictions
    │   │   ├── predict_model.py
    │   │   └── train_model.py
    │   │
    │   └── visualization  <- Scripts to create exploratory and results oriented visualizations
    │       └── visualize.py
    │
    └── tox.ini            <- tox file with settings for running tox; see tox.testrun.org

--------

<p><small>Project based on the <a target="_blank" href="https://drivendata.github.io/cookiecutter-data-science/">cookiecutter data science project template</a>. #cookiecutterdatascience</small></p>

------

Project Flow
------------

**<u>Note:</u>** All descriptions of text files mentioned here are in a README.md in `data/interim` where the text files themselves are located. A lot of text files that are output aren't used in subsequent codes but were created for tracking/visualization purposes. All scripts, unless otherwise stated, can be found in `src/data/main`.

1. prelim_checks_1.py

   | Input       | Original ECO catalog in `/data/raw`      |
   | ----------- | ---------------------------------------- |
   | **Process** | **Adds a filter column, carries out an exposure time check and also whether or not there is imaging in at least 2 filters** |
   | **Output**  | **`Obj_arr.txt`, creates a `URLs` folder where URLs for all images per ECOID are written and kept as text files** |

   This is the first script that was run using just the ECO catalog. Exposure times were calculated using [ETC](http://etc.stsci.edu/etc) per filter for each instrument that was used. The constraint was that the exposure time added up over all images in a particular filter had to meet/exceed the exposure time limit calculated. 208 objects passed.

   ​

2. prelim_checks_2.py

   | Input       | `src/data/mods_prelim_checks_2`, original ECO catalog, `Obj_arr.txt` from `prelim_checks_1.py` (PC1) and `pget_ECO.py` |
   | ----------- | ---------------------------------------- |
   | **Process** | **Takes images that passed PC1 and checks if the galaxy is in the image. After this check has passed the exposure time check and 2-filter check are both repeated in case some images are lost if the object was not in them** |
   | **Output**  | **`obj_rev.txt` is created per object and is a list of all the images that passed this second round of preliminary checks, `goodObj.txt`, `badObj.txt` and `Error_prelim2.txt` ** |

   This is the second script that was run. `pget_ECO.py` is the script that is used to download the images using the URLs in the folders created by `prelim_checks_1.py`. 110 objects passed and 98 failed.

   ​

3. ECO_format.py

   | Input       | Original ECO catalog in `data/raw`       |
   | ----------- | ---------------------------------------- |
   | **Process** | **Adds filter column and new_filename column (without 'http…=' part) permanently and saves these changes to a new file that is used as the ECO catalog from here on** |
   | **Output**  | **`ECO_formatted.txt`**                  |

   This is the third (short) script that reformats the original ECO catalog.

   ​

4. prelim_checks_3.py

   | Input       | `goodObj.txt` from `prelim_checks_2.py` (PC2) and `ECO_formatted.txt` from `ECO_format.py` |
   | ----------- | ---------------------------------------- |
   | **Process** | **Same as PC2 except there is an additional check that is required of whether or not the pixel value at the RA and DEC of galaxy is 0 or not. Exposure time and 2-filter checks are repeated since some images were lost. ** |
   | **Output**  | **`Error_prelim3.txt`, `Prelim3_results.txt`, `Prelim3_results_good.txt`, `Prelim3_results_bad.txt`, `goodObjv2.txt`, and `good_images.txt` per object of all the images that passed this round of checks.** |

   This check was carried out because there were some images where even though the object was in the frame, part of the image was rotated as shown below on the right. The green circle shows the pixel position of the galaxy in question given its RA and DEC and you can see that the object is in the frame but in an area where pixel values are 0. Hence, this check. 

   ![Picture1](/Users/asadm2/Desktop/Picture1.png)

   You may wish to combine all these checks into one script but that will require handling where the text files are read between all codes.

   ​

5. swarp.py

   | Input       | `src/data/mods_swarp.py`, `ECO_formatted.txt`, `goodObjv2.txt` from `prelim_checks_3.py` (PC3) and images per filter per object. |
   | ----------- | ---------------------------------------- |
   | **Process** | **Per filter per object: image units are checked and zeropoints are rescaled (30) if image isn't in counts/s. Headers and data are updated and written to new versions of the images. Swarp is used to combine images per filter into coadds. Finally, the percentage of a coadd that is blank is calculated.** |
   | **Output**  | **New versions of old fits images (`…_test_…`), `coadd.fits`, `coadd.weight.fits`, `swarp.xml`, `percent_blankv2.txt`, `Error_swarp.txt` and `img_list_testfil.txt` per object which is a text file of all the image names from new versions that are created. ** |

   Originally, the blank percentage was calculated by seeing how many pixels had values that were less than/equal to 0 but then we changed it to only those values that are equal to 0. Hence, the second version of `percent_blank.txt`. The `Error_swarp.txt` file keeps a track of all errors that are raised while running this script and should be checked every time swarp is run. [(SWarp manual)](https://www.astromatic.net/pubsvn/software/swarp/trunk/doc/swarp.pdf)

   ------

   **<u>IMPORTANT!</u>** Any time swarp is re-run all fits files that will be re-written will throw an error that is written to `Error_swarp.txt` which is that 'the file already exists.' So, it is important to delete all fits images between runs of swarp.py. You will probably want to write a short script for this. 

   ------

6. coadd_multiply.py

   | Input       | Text file of 'good' objects and `percent_blank.txt` from `swarp.py` |
   | ----------- | ---------------------------------------- |
   | **Process** | **Keeps track of 100% blank coadds and checks, if after removing these coadds, that there are coadds in at least 2 filters. If not, the object can no longer be used otherwise the data of all remaining coadds is multiplied and blank percentage is calculated again but on this combined coadd.** |
   | **Output**  | **`percent_blank_combinedcoaddv3.txt`**  |

   Using `percent_blank.txt` to keep a track of how many coadds were 100% blank, we wanted to see how many of our combined coadds would be blank. Multiplying was the only way a pixel value of 0 would propagate when combining images. The objects that didn't pass this check are manually noted and used in the following script as `[bad_Obj]`. This script did NOT modify the fits images and data in any way. It was merely a check. 

   ​

7. coadd_add.py

   | Input       | Same as `coadd_multiply.py` and manually added `bad_Obj` array from `coadd_multiply.py` |
   | ----------- | ---------------------------------------- |
   | **Process** | **Same as `coadd_multiply.py` except coadds are <u>added</u> together** |
   | **Output**  | **`comb_coadd.fits` per object**         |

   This script is used to combine all coadds per object. This combined coadd is used as the detection image in SExtractor's dual mode. 

   ​

8. sextractor.py

   | Input       | `comacandidates.txt` (manually created from `/src/visualization/eco_sample.py`), `ECO_formatted.txt`, `ECOphot_final.txt`, `comb_coadd.fits`, `acs_f475w` & `acs_f814w` `coadds.fits` per object. |
   | ----------- | ---------------------------------------- |
   | **Process** | **Run SExtractor on subset of Coma that has imaging in f475w and f814w. Petrosian magnitudes are obtained, zeropoints are calculated and applied to petro mags. Lupton transformations are used to compare SE mags to closest equivalent SDSS band magnitude in ECO photometry catalog. (r and i band) 2 different methods of matching between the object to its SE detection are used and tested. (More info in script) ** |
   | **Output**  | **Depending on which method was used: `magnitude_errors_degmatch.txt`, `magnitude_errors_pixelmatch.txt`, `catmatch_separation.txt`. Figures: `/reports/figures/calcrmag_catrmag_coma.png`, `/reports/figures/calcimag_catimag_coma.png` and `reports/figures/catmatch_separation.png`. ** |

   The original goal was to be able to compare our magnitudes to the [Hammer et. al 2010](https://arxiv.org/pdf/1005.3300.pdf) catalog of the Coma cluster as a test of our photometry. The magnitudes weren't matching post swarp.py which is why rescaling the zeropoint was necessary. As commented in the script, the calculations of zeropoints and their application to the petro mags are now outdated since the rescaling was applied. Now, you just have to add the new zeropoint directly to the magnitudes that SE returns. 

   [Lupton transformations](http://www.sdss3.org/dr8/algorithms/sdssUBVRITransform.php#Lupton2005)

   [SExtractor manual](http://astroa.physics.metu.edu.tr/MANUALS/sextractor/Guide2source_extractor.pdf)

   ------

   We then decided to switch gears and try and get an estimate of how many GCs we can detect in a subsample. All the following scripts are in `src/visualization`.

   ​

9. eco_sample.py

   | Input       | `eco_dr1.txt`, `ECOphot_final.txt` and `goodObjv2.txt` |
   | ----------- | ---------------------------------------- |
   | **Process** | **Plot current sample (Figure 1) and use equation 1 in Zaritsky et al. 2015 to get an estimate of number of globular clusters (GCs) ** |
   | **Output**  | **N/A**                                  |

   [Zaritsky et al. 2015](https://arxiv.org/pdf/1511.05608.pdf)

10. MC.py

    | Input       | Same as `eco_sample.py`                  |
    | ----------- | ---------------------------------------- |
    | **Process** | **Run 1000 realizations of MC using gaussian distributions for number of GCs given stellar mass of a galaxy and gaussian distribution for magnitudes to get total luminosity of the GC system. Then using a M/L ratio, get a mass estimate of the GC system and ultimately a halo mass. This was tested for <u>ECO12028</u> which is the smallest halo in our `single_halos` subsample from `eco_sample.py`.** |
    | **Output**  | **N/A**                                  |

    The sources for mean and sigma of both distributions are mentioned in script's comments. This script has just been run for the lowest mass halo which is the last thing that was done. The luminosity function for GCs (Figure 2) and the dependence of sigma on b-band magnitude (Equation 18) both come from [Jordán et al. 2007](https://arxiv.org/pdf/astro-ph/0702496.pdf)