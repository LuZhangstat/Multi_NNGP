


Structure
---------
|Folder Name |Subfolder Name | File Name|     Intro            |
|:----- |:----- |:----- |:-------------- |
|**data**|***rawdata***| readdata.R | download rawdata|
|    |***EDA***| EDA.R | exploratory data analysis|                               |
|**projects**|***Conj***| data_conj.ipynb | Precalculation for conjugate models|
|        |      | Bayesian_linear_real.ipynb | Bayesian linear model |
|        |      | Multi_conj_latent_real_parallel.ipynb | Multivariate conjugate latent model |
|        |      | Multi_conj_latent_real_summary.ipynb | Summary for multivariate conjugate latent model |
|        |      | Multi_conj_resp_real_parallel.ipynb | Multivariate conjugate response model |
|        |      | Multi_conj_resp_real_summary.ipynb | Summary for multivariate conjugate response model |
|        |***misalign_response***|Conj_res_misalign_data.ipynb | Precalculation for conjugate response model with misalignment|
|        |      |Bayesian_linear_real_small.ipynb | Bayesian linear model for subset of whole data|
|        |      | *misalign_real.ipynb | Multivariate conjugate response model with misalignment |
|        |      | *misalign_summary.ipynb | Summary for onjugate response model with misalignment |
|        |***BSLMC***| data_BSLMC.ipynb | Precalculation for BSLMC|
|        |      | BSLMC_forloop*.ipynb | BSLMC model |
|        |      | BSLMC_summary*.ipynb | Summary for BSLMC model |
|        |***Factor_BSLMC***|data_Factor_BSLMC.ipynb | Precalculation for factor BSLMC with diagonal Sigma|
|        |      | BSLMC_Factor_real-RAM*.ipynb | factor BSLMC model with diagonal Sigma |
|        |      | BSLMC_Factor_real-summary*.ipynb | Summary for factor BSLMC model with diagonal Sigma |


Notes
---------
* We use "*" as a shorthand of file name in the above table 
* One need to compile the library in folder "../julia-R-nn-ccall2" before runing codes in this folder
* One need to creat folders for data, results and pictures based on the code.


