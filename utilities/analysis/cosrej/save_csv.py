import ROOT
import pandas as pd
import numpy as np

t_signal = ROOT.RDataFrame("doms", "bdtinput_signal_aframe.root")
t_bkg = ROOT.RDataFrame("doms", "bdtinput_bkg_aframe.root")

df_signal = pd.DataFrame(t_signal.AsNumpy())
df_bkg = pd.DataFrame(t_bkg.AsNumpy())

colagg = {}
info_cols = ['rock_wgt', 'vtxX', 'vtxY', 'vtxZ']
for col in df_signal.columns:
    if col == "event_id":
        continue
    if col not in info_cols:
        colagg[col] = list
    else:
        colagg[col] = np.mean

df_signal = df_signal.groupby("event_id").agg(colagg).reset_index()
df_bkg = df_bkg.groupby("event_id").agg(colagg).reset_index()

max_len = max(df_signal["dom_x"].apply(len).max(), df_bkg["dom_x"].apply(len).max())

def pad_trunc(l, length, pad_value=0):
    if len(l) >= length:
        return l[:length]
    return l + [pad_value] * (length - len(l))

# Apply to all columns except event_id and various truth info columns
for col in df_signal.columns:
    if col == "event_id" or col in info_cols:
        continue
    df_signal[col] = df_signal[col].apply(lambda l: pad_trunc(l, max_len, pad_value=-1))
    df_bkg[col] = df_bkg[col].apply(lambda l: pad_trunc(l, max_len, pad_value=-1))


def df_to_numpy(df):
    arrays = []
    for c in df.columns:
        if c == 'event_id' or c in info_cols:
            continue
        if type(df[c][0]) is not list:
            arrays.append(np.array(df[c].to_list()).reshape(-1, 1))
        else:
            arrays.append(np.array(df[c].to_list()))
    return np.hstack(arrays)

X_signal = df_to_numpy(df_signal)  # shape (n_signal,  L * n_cols)
X_bkg = df_to_numpy(df_bkg)  # shape (n_bkg,  L * n_cols)
df_signal_info = df_signal[info_cols]
df_bkg_info = df_bkg[info_cols]

print(X_signal.shape, X_bkg.shape)

np.save('signal_aframe.npy', X_signal)
np.save('bkg_aframe.npy', X_bkg)
df_signal_info.to_csv('signal_aframe_info.csv', index=False)
df_bkg_info.to_csv('bkg_aframe_info.csv', index=False)

with open('info_aframe.txt', 'w') as info:
    info.write('Columns : %s'%str(list(df_signal.columns)))
    info.write('\nArray Length: %d'%max_len)
