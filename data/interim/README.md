### README file for all txt files

**<u>Note:</u>** Some text files are only stored locally and may not be present in the GitHub repo.

| File                          | Description                              |
| ----------------------------- | ---------------------------------------- |
| badObj                        | List of objects that didn't pass exposure time and multiple filter checks |
| ECO_formatted                 | Original ECO catalog formatted to have a 'filters' and 'new_filename' column |
| goodandbadObj                 | Both goodObj and badObj combined to find the missing objects during exposure time and multiple filter checks |
| goodObj                       | List of good objects after running prelim_checks_2.py |
| Obj_arr                       | Original list of good objects after running prelim_checks_1.py |
| percent_emptyfil              | List of coadded images per object and per filter and what percent empty they are i.e. pixels have either 0 or negative values |
| Error_swarp                   | Any errors when running swarp per filter per object |
| percent_blank_combinedcoadd   | After adding all coadds per filter and taking into account all pixel values less than OR exactly 0 and counting them as 'blank' |
| percent_blank_combinedcoaddv2 | After adding all coadds per filter and taking into account all pixel values that are ONLY exactly 0 and counting them as 'blank' |
| percent_blank_combinedcoaddv3 | Ran coadd_multiply.py again              |
| percent_blank                 | Calculating per object per filter and taking into account all pixel values less than OR exactly 0 and counting them as 'blank' |
| percent_blankv2               | Calculating per object per filter and taking into account all pixel values that are ONLY exactly 0 and counting them as 'blank' |
| Error_prelim3                 | Errors running prelim_checks_3.py        |
| ECOphot_final                 | ECO photometry catalog                   |
| eco_dr1                       | Original ECO catalog from website        |
| Prelim3_results_bad           | Had to repeat this because there were images where the object was there but in the black region because the image was rotated. |
| Prelim3_results_good          | Had to repeat this because there were images where the object was there but in the black region because the image was rotated. |
| goodObjv2                     | Most recent list of good objects after running prelim_checks_3.py |
| Prelim3_results               | Writing results for each object: number of good images, total number of images. Had to repeat this because there were images where the object was there but in the black region because the image was rotated. |
| percent_emptyfil_only0        | List of coadded images per object and per filter and what percent empty they are i.e. pixels have a value of EXACTLY 0 |
| catmatch_separation           | Separation between real object and SE detection when using match_to_catalog method in sextractor.py |
| comacandidates                | Manually created from eco_sample.py to be the ECOIDs of all objects with group halo mass greater than 10^14^ solar masses |
|                               |                                          |









