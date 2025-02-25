import uproot
import matplotlib.pyplot as plt
import numpy as np
import argparse

def main():
    args = GetArguments()
    tfile = uproot.open(args.file)
    outputTree = tfile['output']
    metaTree = tfile['meta']
    # Load meta table
    meta = metaTree.arrays(metaTree.keys(), library="np")
    output = outputTree.arrays(outputTree.keys(), library="np")
    # Plot Time and charge of pmts

    ### Time (ns)
    timeoptions = { "bins": np.arange(-20, 120, 0.2),
                    "histtype": 'step' }
    plt.hist(np.concatenate(output['hitPMTTime']), label="hitPMTTime", **timeoptions)
    plt.hist(np.concatenate(output['hitPMTDigitizedTime'])-args.offset, label="hitPMTDigitizedTime", **timeoptions)
    plt.xlabel("Time (ns)")
    plt.ylabel("Count")
    plt.legend()
    plt.show()

    ## Charge (pC)
    chargeoptions = { "bins": np.arange(-2, 20, 0.1),
                      "histtype": 'step' }
    plt.hist(np.concatenate(output['hitPMTCharge']), label="hitPMTCharge", **chargeoptions)
    plt.hist(np.concatenate(output['hitPMTDigitizedCharge']), label="hitPMTDigitizedCharge", **chargeoptions)
    plt.xlabel("Charge (pC)")
    plt.ylabel("Count")
    plt.legend()
    plt.show()

def GetArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("file")
    parser.add_argument("--offset", default=67.7)
    return parser.parse_args()

if __name__ == '__main__':
    main()
