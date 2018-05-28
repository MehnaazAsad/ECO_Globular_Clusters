### README file for all txt files

| File                          | Description                              |
| ----------------------------- | ---------------------------------------- |
| badObj                        | List of objects that didn't pass exposure time and multiple filter checks |
| ECO_formatted                 | Original ECO catalog formatted to have a 'filters' and 'new_filename' column |
| ECO_idurl                     | ECOIDs with their respective image URLs (one of the outputs from WIP) |
| ECO_url                       | List of all URLs (one of the outputs from WIP) |
| goodandbadObj                 | Both goodObj and badObj combined to find the missing objects during exposure time and multiple filter checks |
| goodObj                       | List of objects that did pass exposure time and multiple filter checks |
| percent_emptyfil              | List of coadded images per object and per filter and what percent empty they are i.e. pixels have either 0 or negative values |
| Error_swarp                   | Any errors when running swarp per filter per object |
| percent_blank_combinedcoadd   | After adding all coadds per filter and taking into account all pixel values less than OR exactly 0 and counting them as 'blank' |
| percent_blank_combinedcoaddv2 | After adding all coadds per filter and taking into account all pixel values that are ONLY exactly 0 and counting them as 'blank' |
| percent_blank_combinedcoaddv3 | Ran coadd_multiply.py again              |
| percent_blank                 | Calculating per object per filter and taking into account all pixel values less than OR exactly 0 and counting them as 'blank' |
| percent_blankv2               | Calculating per object per filter and taking into account all pixel values that are ONLY exactly 0 and counting them as 'blank' |
| Error_prelim3                 | Errors running prelim_checks_3.py               |
| ECOphot_final                 | ECO photometry catalog                   |
| eco_dr1                       | Original ECO catalog from website        |
| Prelim3_results_bad           | Had to repeat this because there were images where the object was there but in the black region because the image was rotated. |
| Prelim3_results_good          | Had to repeat this because there were images where the object was there but in the black region because the image was rotated. |
| goodObjv2                     | New version after getting results above. |
| Prelim3_results               | Writing results for each object: number of good images, total number of images. Had to repeat this because there were images where the object was there but in the black region because the image was rotated. |
| percent_emptyfil_only0        | List of coadded images per object and per filter and what percent empty they are i.e. pixels have a value of EXACTLY 0 |















