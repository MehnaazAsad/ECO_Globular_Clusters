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
| Error_swarpfil                | Any errors when running swarp per filter per object |
| percent_blank_combinedcoadd   | After adding all coadds per filter and taking into account all pixel values less than OR exactly 0 and counting them as 'blank' |
| percent_blank_combinedcoaddv2 | After adding all coadds per filter and taking into account all pixel values that are ONLY exactly 0 and counting them as 'blank' |
| percent_blank                 | Calculating per object per filter and taking into account all pixel values less than OR exactly 0 and counting them as 'blank' |
| percent_blankv2.txt           | Calculating per object per filter and taking into account all pixel values that are ONLY exactly 0 and counting them as 'blank' |
| Error_objinimg.txt            | Errors running ObjinImg.py               |
| Obj_in_Img_results.txt        | Writing results for each object: number of good images, total number of images |
| Expfil2_results_bad.txt       | Writing results of running Obj_in_Img.py of all the bad objects for which there aren't enough filters |
| Expfil2_results_good.txt      | Writing results of running Obj_in_Img.py of all the good objects: number of good images,total number of images,number of filters |
| ECOphot_final.txt             | ECO photometry catalog                   |

















