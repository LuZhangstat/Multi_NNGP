


Structure
---------
|Folder Name |Subfolder Name | File Name|     Intro            |
|:----- |:----- |:----- |:-------------- |
|**data**|***rawdata***| readdata.R | download rawdata|
|    |***EDA***| EDA.R | exploratory data analysis|                               |
|**projects**|***Conj*** (no long presented in the paper)| data_conj.ipynb | Precalculation for conjugate models|
|        |      | Bayesian_linear_real.ipynb | Bayesian linear model |
|        |      | Multi_conj_latent_real_parallel.ipynb | Multivariate conjugate latent model |
|        |      | Multi_conj_latent_real_summary.ipynb | Summary for multivariate conjugate latent model |
|        |      | Multi_conj_resp_real_parallel.ipynb | Multivariate conjugate response model |
|        |      | Multi_conj_resp_real_summary.ipynb | Summary for multivariate conjugate response model |
|        |***misalign_response*** (no long presented in the paper)|Conj_res_misalign_data.ipynb | Precalculation for conjugate response model with misalignment|
|        |      |Bayesian_linear_real_small.ipynb | Bayesian linear model for subset of whole data|
|        |      | *misalign_real.ipynb | Multivariate conjugate response model with misalignment |
|        |      | *misalign_summary.ipynb | Summary for onjugate response model with misalignment |
|        |***BLMC***| data_BSLMC.ipynb | Precalculation for BSLMC|
|        |      | BSLMC_forloop*.ipynb | BLMC model |
|        |      | BSLMC_summary*.ipynb | Summary for BLMC model |
|        |      | linear_forloop.ipynb | Bayesian linear model in Section 5|
|        |***Factor_BLMC***|data_Factor_BSLMC.ipynb | Precalculation for BLMC with diagonal Sigma|
|        |      | BSLMC_Factor_real-RAM*.ipynb | BLMC model with diagonal Sigma |
|        |      | BSLMC_Factor_real-summary*.ipynb | Summary for BLMC model with diagonal Sigma |


Notes
---------
* We use "*" as a shorthand of file name in the above table 
* One need to compile the library in folder "../julia-R-nn-ccall2" before runing codes in this folder


