import os
import shutil
import glob

KEEP_EVERY = 100
KEEP_LAST = 2
WORKDIR = "./NMC047_deposition"

DRY_RUN = False

file_types_to_remove = [
"final-coordinates",
"checkpoint",
"deposition-log",
"topology",
"input-coordinates",
"control",
"restraints",
"stdout",
"stderr",
"tpr",
"log",
"trajectory",
"energy",
]

def clean_workdir():
    run_ids = [int(f.split("_")[-1].split(".")[0]) for f in os.listdir(os.path.join(WORKDIR, "final-coordinates")) if len(f.split("_"))>2]
    to_keep = []
    to_remove = []
    last_run = max(run_ids)
    for i in sorted(run_ids):
        if not i % KEEP_EVERY or i > last_run-KEEP_LAST:
            to_keep.append(i)
        else:
            to_remove.append(i)
    print("runs to keep: " + str(to_keep))
    print("runs to remove:" + str(to_remove))
    if not DRY_RUN:
        print("This is NOT a dry run, files WILL be deleted.")
    for i in to_remove:
        print("cleaning run {}".format(i))
        for file_type in file_types_to_remove:
            for f in glob.glob(os.path.join(WORKDIR, file_type, "*_{}_{}.*".format(file_type, i))):
                if DRY_RUN:
                    print("\twould be deleting: {}".format(f))
                else:
                    print("\tdeleting: {}".format(f))
                    os.remove(f)

if __name__=="__main__":
    clean_workdir()
