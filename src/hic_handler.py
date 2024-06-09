def line_counter(file: str):
    c = 0
    with open(file, "r") as fd:
        for _ in fd:
            c += 1
        fd.close()
    return c

def load_contracts(bed_file: str, matrix_file: str):
    num_entry = line_counter(bed_file)
    bed_arr = [None for _ in range(num_entry + 1)] # 1-based
    name_idx_dict = {}

    with open(bed_file, "r") as bed_fd:
        for line in bed_fd:
            name, start, end, idx = line.strip().split("\t")
            bed_arr[int(idx)] = (name, int(start), int(end))
            if name not in name_idx_dict:
                name_idx_dict[name] = []
            name_idx_dict[name].append(int(idx))
        bed_fd.close()

    name_st_dict = {}
    for name, arr in name_idx_dict.items():
        s, t = arr[0], arr[-1]
        name_st_dict[name] = (s, t)
    inter_contacts = {}
    intra_contacts = {}
    for name in name_st_dict.keys():
        intra_contacts[name] = []
    with open(matrix_file, "r") as mat_fd:
        for line in mat_fd:
            uidx, vidx, val = line.strip().split("\t")
            link = (int(uidx), int(vidx), int(round(float(val))))
            uname = bed_arr[int(uidx)][0]
            vname = bed_arr[int(vidx)][0]
            if uname == vname:
                # intra contacts
                intra_contacts[uname].append(link)
            else:
                # inter contacts
                if (uname, vname) not in inter_contacts:
                    inter_contacts[(uname, vname)] = []
                inter_contacts[(uname, vname)].append(link)
        mat_fd.close()
    return bed_arr, name_st_dict, inter_contacts, intra_contacts
