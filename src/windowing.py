import sys
import numpy
import matplotlib.pyplot as plt
import seaborn as sns

def inspect_window(name_st_dict: dict, inter_contacts: dict, hidx: dict, cname: str, bs: str, windex: str):
    window_dict = {}
    for (uname, vname), links in inter_contacts.items():
        if uname not in hidx or vname not in hidx:
            continue
        if uname != cname and vname != cname:
            continue
        if uname == cname:
            if vname not in window_dict:
                window_dict[vname] = []
            for (uidx, vidx, val) in links:
                ruidx = uidx - bs + 1
                if ruidx == windex:
                    window_dict[vname].append((vidx-name_st_dict[vname][0]+1, val))
        else:
            if uname not in window_dict:
                window_dict[uname] = []
            for (uidx, vidx, val) in links:
                rvidx = vidx - bs + 1
                if rvidx == windex:
                    window_dict[uname].append((uidx-name_st_dict[uname][0]+1, val))
    print(sum(sum(link[1] for link in links) for links in window_dict.values()))
    return

def plot_window(name_st_dict: dict, inter_contacts: dict, hidx: dict, out_dir: str, window_size: str, visual=False):
    for cname, (bs, bt) in name_st_dict.items():
        if cname not in hidx:
            continue
        print("==============================")
        print(f"accessing {cname}, absolute window index: [{bs}...{bt}]")
        cwindow = dict.fromkeys(range(1, bt - bs + 2), 0)
        for (uname, vname), links in inter_contacts.items():
            if uname not in hidx or vname not in hidx:
                continue
            if uname != cname and vname != cname:
                continue
            if uname == cname:
                for (uidx, vidx, val) in links:
                    ruidx = uidx - bs + 1
                    cwindow[ruidx] += val
            else:
                for (uidx, vidx, val) in links:
                    rvidx = vidx - bs + 1
                    cwindow[rvidx] += val

        max_tick = max(cwindow.keys(), key=lambda x: cwindow[x])
        print(cwindow.items())
        print(max_tick, cwindow[max_tick])
        if cname == "chr12B":
            inspect_window(name_st_dict, inter_contacts, hidx, cname, bs, max_tick)
        if visual:
            _, ax = plt.subplots(1, 1, figsize=(30, 15))
            ind = numpy.arange(1, len(cwindow)+1)
            palette = sns.color_palette("husl", len(cwindow))
            plt.bar(ind, list(cwindow.values()), color=palette)
            ticks = [""]*len(cwindow.keys())
            ticks[max_tick - 1] = f"{max_tick}*{window_size}bp"
            plt.xticks(ind, ticks)
            plt.savefig(f"{out_dir}/{cname}.png")
            plt.close()

    print("==============================")