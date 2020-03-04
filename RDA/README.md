


Roadmap
---------
|Folder Name |     Intro            |
|:------ |:----------- |
|sim1-3| old simulations, not presented in paper|
|sim4| simulation example 1 |
|sim_Factor| simulation example 2 |


Structure
------------------
+-- data<br />
|   &nbsp;&nbsp;       +-- (step 1)rawdata: run "readdata.R" to download rawdata <br />
|   &nbsp;&nbsp;       +-- (step 2)EDA: run "EDA.R" to do exploratory data analysis <br />
|    &nbsp; &nbsp;  &nbsp;   &nbsp;         
|<br />
+-- projects: R code files for full GP, response NNGP and latent NNGP.<br />
|   &nbsp;&nbsp;       +-- Conj: run "readdata.R" to download rawdata <br />
|   &nbsp;&nbsp; &nbsp;&nbsp; +-- (step1)
|   &nbsp;&nbsp;       +-- BSLMC: run "EDA.R" to do exploartory data analysis <br />
|   &nbsp;&nbsp;       +-- misalign_response: run "readdata.R" to download rawdata <br />
|   &nbsp;&nbsp;       +-- Factor_BSLMC: run "EDA.R" to do exploartory data analysis <br />
|<br />


Notes
---------
* One need to compile the library in folder "./julia-R-nn-ccall2" before runing codes in this folder




