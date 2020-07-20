"""
This file should contain all the references needed to read from a file and store
the data properly. For dictionaries containing figure labels and such see
the Visualization module.
"""


# Known raw data dictionary
raw_data_dict = dict()
raw_data_dict["T0"] = ["#T0(K)"]
raw_data_dict["I"] = ["I(A)"]
raw_data_dict["R+_0"] = ["R+_0(V)"]
raw_data_dict["R+_Q"] = ["R+_Q(V)"]
raw_data_dict["R-_0"] = ["R-_0(V)"]
raw_data_dict["R-_Q"] = ["R-_Q(V)"]
raw_data_dict["dTy_0"] = ["dTy_0(V)"]
raw_data_dict["dTy_Q"] = ["dTy_Q(V)"]
raw_data_dict["dTabs_0"] = ["Tabs_0(V)"]
raw_data_dict["dTabs_Q"] = ["Tabs_Q(V)"]
raw_data_dict["dTx_0"] = ["dTx_0(V)"]
raw_data_dict["dTx_Q"] = ["dTx_Q(V)"]

# Known measurements dictionary
measurements_dict = dict()
measurements_dict["T0"] = ["T0(K)", "T0 (K)"]
measurements_dict["T_av"] = ["T_av(K)", "Taverage(K)", "T (K)"]
measurements_dict["Tp"] = ["T+(K)", "T+ (K)"]
measurements_dict["Tm"] = ["T-(K)", "T- (K)"]
measurements_dict["dTx"] = ["dTx(K)", "dTx (K)"]
measurements_dict["kxx"] = ["kxx(W/Km)", "k_xx(W/Km)", "Kxx (W / K m)"]
measurements_dict["kxy"] = ["kxy(W/Km)", "k_xy(W/Km)", "Kxy (W / K m)"]
measurements_dict["kxy/kxx"] = ["kxy/kxx(%)"]
measurements_dict["dTy/dTx"] = ["dTy/dTx(%)"]
measurements_dict["dTy"] = ["dTy(K)", "dTy (K)"]
measurements_dict["I"] = ["I(A)", "I (A)"]
measurements_dict["dTabs"] = ["dTabs", "dT_abs"]
measurements_dict["kxx/T"] = ["kxx/T"]
measurements_dict["Resistance"] = ["Resistance"]
measurements_dict["dTx/T"] = ["dTx/T"]
measurements_dict["Tp_Tm"] = ["Tp_Tm"]
measurements_dict["I_fit"] = ["I_fit"]
measurements_dict["T0_fit"] = ["T0_fit"]

# Known parameters dictionary
parameters_dict = dict()
parameters_dict["H"] = ["H"]
parameters_dict["w"] = ["w"]
parameters_dict["t"] = ["t"]
parameters_dict["L"] = ["L"]
parameters_dict["mount"] = ["mount"]
parameters_dict["sample"] = ["Sample", "sample"]
parameters_dict["date"] = ["Date", "date"]

