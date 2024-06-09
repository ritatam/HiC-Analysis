import os
import sys

import numpy

from hic_handler import load_contracts
from view_circos import draw_circos_plot
from windowing import plot_window


if __name__ == "__main__":
    if len(sys.argv) != 8:
        print("Generate HiC statistics and circos plot")
        print(f"{sys.argv[0]} <hap1.lst> <hap2.lst> <na.lst> <bed> <matrix> <outdir> <w>")
        print("<hap1.lst> path to a file consists the list of contig id prefixed with > from haplotype 1, one per line")
        print("<hap2.lst> path to a file consists the list of contig id prefixed with > from haplotype 2, one per line")
        print("<na.lst>   path to a file consists the list of contig id prefixed with > from unclassified, one per line")
        print("<bed>      path to HiC-Pro bed file")
        print("<matrix>   path to HiC-Pro matrix file")
        print("<outdir>   output directory")
        print("<w>        HiC-Pro window size")
        
        sys.exit(1)
    
    [_, hap1file, hap2file, na_file, bed_file, mat_file, out_dir] = sys.argv[:7]

    window_size = int(sys.argv[7])

    if out_dir.endswith('/'):
        out_dir = out_dir[:-1]
    
    os.makedirs(out_dir, exist_ok=True)
    out_wdir = f"{out_dir}/window_plot/"
    os.makedirs(out_wdir, exist_ok=True)

    # parse input
    hap_dict = {}
    hap_rec = {1: [], 2: [], 3: []}
    hidx = {}
    header = []
    hcount = 0
    for (i, hap_file) in [(1, hap1file), (2, hap2file), (3, na_file)]:
        with open(hap_file, 'r') as fd:
            for line in fd:
                hap = line.strip()[1:]
                assert hap not in hap_dict
                hap_dict[hap] = i
                hap_rec[i].append(hap)
                header.append(hap)
                hidx[hap] = hcount
                hcount += 1
            fd.close()

    # load HiC matrix
    _, name_st_dict, inter_contacts, intra_contacts = load_contracts(bed_file, mat_file)
    hic_mat = numpy.zeros((hcount, hcount), dtype=numpy.int64)

    total_link = 0

    hap1_link = 0
    hap2_link = 0
    ctg1_internal = 0
    ctg2_internal = 0

    crossing_links = 0

    ext_link = 0
    ext_internal = 0

    # unclassified to hap 1/2
    ext1_link = 0 
    ext2_link = 0

    # only consider the classified ones
    ## within chromosome
    cis_link = 0
    ## across chromosome
    trans_link = 0

    for uname, links in intra_contacts.items():
        if uname not in hidx:
            continue
        sum_weight = sum(val for (_, _, val) in links)
        hic_mat[hidx[uname], hidx[uname]] = sum_weight
        total_link += sum_weight
        cis_link += sum_weight
        if hap_dict[uname] == 1:
            ctg1_internal += sum_weight
        elif hap_dict[uname] == 2:
            ctg2_internal += sum_weight
        else:
            ext_internal += sum_weight

    for (uname, vname), links in inter_contacts.items():
        if uname not in hidx or vname not in hidx:
            continue
        sum_weight = sum(val for (_, _, val) in links)
        hic_mat[hidx[uname], hidx[vname]] += sum_weight
        hic_mat[hidx[vname], hidx[uname]] += sum_weight
        total_link += sum_weight

        if uname[:-1] == vname[:-1]:
            cis_link += sum_weight
        else:
            trans_link += sum_weight

        if hap_dict[uname] != hap_dict[vname]:
            if hap_dict[uname] in [1,2] and hap_dict[vname] in [1,2]:
                crossing_links += sum_weight
            else:
                if hap_dict[uname] == 1 or hap_dict[vname] == 1:
                    ext1_link += sum_weight
                if hap_dict[uname] == 2 or hap_dict[vname] == 2:
                    ext2_link += sum_weight
        else:
            if hap_dict[uname] == 1:
                hap1_link += sum_weight
            elif hap_dict[uname] == 2:
                hap2_link += sum_weight
            else:
                ext_link += sum_weight

    # write matrix to file
    out_hic_mat = f"{out_dir}/hic_mat.csv"
    hic_header = hap_rec[1] + hap_rec[2] + hap_rec[3]
    with open(out_hic_mat, 'w') as fd:
        fd.write(',' + ",".join(hic_header) + '\n')
        for ridx in range(hcount):
            us = hcount + 1
            disp = [header[ridx]]
            for cidx in range(hcount):
                disp.append(str(hic_mat[ridx, cidx]))
            fd.write(",".join(disp) + '\n')            
        fd.close()
    numpy.savetxt(f"{out_dir}/hic_mat.raw.csv", hic_mat, delimiter=",")
    with open(f"{out_dir}/hic_mat_header.csv", 'w') as fd:
        fd.write(",".join(header) + '\n')
        fd.close()

    print(f"Total Link\t{total_link}")
    print(f"CIS Link\t{cis_link}\t{round(cis_link/total_link, 2)*100}%")
    print(f"Trans Link\t{trans_link}\t{round(trans_link/total_link, 2)*100}%")
    print(f"HAP1_Internal\t{ctg1_internal}\t{round(ctg1_internal/total_link, 2)*100}%")
    print(f"HAP2_Internal\t{ctg2_internal}\t{round(ctg2_internal/total_link, 2)*100}%")
    print(f"HAP1 Link\t{hap1_link}\t{round(hap1_link/total_link, 2)*100}%")
    print(f"HAP2 Link\t{hap2_link}\t{round(hap2_link/total_link, 2)*100}%")
    print(f"CROSS HAP Link\t{crossing_links}\t{round(crossing_links/total_link, 2)*100}%")
    print(f"EXT-EXT Link\t{ext_link}\t{round(ext_link/total_link, 2)*100}%")
    print(f"EXT-HAP1 Link\t{ext1_link}\t{round(ext1_link/total_link, 2)*100}%")
    print(f"EXT-HAP2 Link\t{ext2_link}\t{round(ext2_link/total_link, 2)*100}%")
    print(f"EXT_Internal\t{ext_internal}\t{round(ext_internal/total_link, 2)*100}%")


    clen_dict = {}
    for uname, (us, ut) in name_st_dict.items():
        clen_dict[uname] = window_size * (ut-us+1)

    if not os.path.exists(f"{out_dir}/circos_plot.png"):
        draw_circos_plot(clen_dict, hap_dict, hap_rec,
                        window_size, total_link, inter_contacts, 
                        name_st_dict, out_dir)
    
    plot_window(name_st_dict, inter_contacts, hidx, out_wdir, window_size, True)
    
    print(f"output matrix:                {out_dir}/hic_mat.csv")
    print(f"output raw hic matrix:        {out_dir}/hic_mat.raw.csv")
    print(f"output raw hic matrix header: {out_dir}/hic_mat_header.csv")
    print(f"output circos plot:           {out_dir}/circos_plot.png")
    print(f"output window plots:          {out_wdir}/*.png")
    print("Done")

    sys.exit(0)
    
