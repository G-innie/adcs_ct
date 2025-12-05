### How to use
this repo contains mainly wrapper files `gain_search.py` that can be used to look for good gains for the controller for mist. the range of considered $k_p$ is to be initialised in the following lines. The exp part is about how high of an exponential one wnats to consider, whereas the mul part is about how fine you want the gain search to be. So, for example, with $k_p=0.25 \cdot \text{mul}$, and `mul in range(4, 40)`, we consider values such that $1\cdot10^\text{exp}, 1.25\cdot10^\text{exp}, 1.5\cdot10^\text{exp}, ..., 2\cdot10^\text{exp}, ... 9.75\cdot10^\text{exp}$. The $k_d$ has been chosen to be ten times as large as $k_p$, that is $k_d = 10\cdot k_p$. Once a good $k_p$ has been found one can repeat this process to look for a better $k_d$ (but on a smaller scale ofc).

```
kp_vals = []
for exp in range(1, 40):  # 10^4 to 10^40
    for mul in range(4, 40):  # 1*10^exp, 1.25*10^exp, ..., 2*10^exp, ... 9*10^exp
        kp_vals.append(0.25*mul * 10**exp)
kp_vals = np.array(kp_vals)
``` 

Then, one needs to adjust the path or where the files live for every `with open`. Depending on your parameters settings, this can take a while to run, for example for simulating 86400 seconds, one can expect around 20 (real life) seconds per simulation. After the gain search is complete, you can find the results for every $k_p$ value explored in `gain_search_results.txt`, where the file saves in each row the value of $k_p$, the value of $k_d$, the mean, the median, the min and max values of all erros, and the final value is the percentage of errors in the sunlight that have value 20 or less.

Then, we plot these statistics for ease of use. (maybe)

We are looking for gains such that the last value is as close to $100\%$ as possible. It could also be good to pay attention to the max error, since we don't want this error to be very large, even in the eclipse. This process should narrow down the gains of interest, and then individual simulations can be rerun to see which of the remaining gains gives the most hopeful results.  

### gain search vs gain search 2
`gain_search2.py` is the more fancy version of `gain_search.py`. In comparison, it calculates the percentage of the errors in the sunlight that are 20 degrees or less, and the function `short_sci` also can handle higher powers (if i remember right this is mostly for printing and saving but still kinda important).

### gain_serch_results_xxx.txt data description
- 1-3 - made with param values as in this repo, not good values to use
- h1 - made with param values as the default in AOCS_oDyn_sim, good values to use
- 20251204 - also made with param values as the default in AOCS_oDyn_sim but also the main.c has been played around with in the following way: noise bool in line 92 is on, albedo in line 206 is on, actuation in line 214 is 100%, activation flag in line 159 is on.
- 20251205 - same as 20251204 but the albedo threshold has been changed from constant to be the largest element in the last photodiodes as described in fig. 26 in M160-032


### other files
Other files are just a bunch of results saved up, I also saved up plots of some more interesting results in case I want to get back to it later. The param files are here for my convenience when moving between laptops.

### to do:
1. modify the `gain_search_ekf_params.py` to do a search for good values of ekf params.