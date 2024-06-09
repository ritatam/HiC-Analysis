import pycircos

def draw_circos_plot(clen_dict: dict, hap_dict: dict, hap_rec: dict, 
                     window_size: int, _: int,
                    inter_contacts: dict, name_st_dict: dict, out_dir: str):    
    colors = {1: "Red", 2: "Blue", 3: "Green"}
    one_tick = 1000000

    rename = True
    label_mapper = {}
    
    Garc = pycircos.Garc
    Gcircle = pycircos.Gcircle
    circle = Gcircle(figsize=(64, 64))
    for hidx, hids in hap_rec.items():
        for i, hid in enumerate(hids):
            if clen_dict[hid] < one_tick / 2:
                continue
            if rename:
                label_mapper[hid] = f'{i+1}' + {1: 'A', 2: 'B', 3: 'H'}[hidx]
            else:
                label_mapper[uname] = uname
            arc = Garc(arc_id=label_mapper[hid], 
                       size=clen_dict[hid],
                       interspace=0.9,
                       raxis_range=(835, 885),
                       labelposition=60,
                       labelsize=48,
                       label_visible=True,
                       facecolor=colors[hidx])
            circle.add_garc(arc)

    circle.set_garcs()

    for arc_id in circle.garc_dict:
        circle.tickplot(
            arc_id, raxis_range=(885, 900), tickinterval=one_tick, ticklabels=None
        )
    
    chord_records = {"Red": [], "Blue": [], "Green": []}
    total_crossing_links = 0
    for (uname, vname), links in inter_contacts.items():
        if uname not in label_mapper or vname not in label_mapper:
            continue

        color = ""
        if hap_dict[uname] == hap_dict[vname]:
            color = colors[hap_dict[uname]]
        else:
            color = "Green"

        us, _ = name_st_dict[uname]
        vs, _ = name_st_dict[vname]
        for (uidx, vidx, val) in links:
            ruidx = uidx - us + 1
            rvidx = vidx - vs + 1

            src = (label_mapper[uname], ruidx*window_size, (ruidx+1)*window_size, 835)
            dest = (label_mapper[vname], rvidx*window_size, (rvidx+1)*window_size, 835)

            chord_records[color].append((src, dest, val))
            total_crossing_links += val

    for color in ["Red", "Blue", "Green"]:
        alpha = float(sum(x for (_,_,x) in chord_records[color]) / total_crossing_links)
        print(color, alpha)
        for (src, dest, val) in chord_records[color]:
            circle.chord_plot(src, dest, facecolor=(color, alpha))
        
    circle.figure.savefig(f"{out_dir}/circos_plot.png", dpi=200.0)

    if rename:
        for oid, nid in label_mapper.items():
            print(nid, oid)