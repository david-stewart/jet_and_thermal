import ROOT as rt
import os
import pandas as pd
from glob import glob as glob
import numpy as np
import json
import sys


# def make_single_bin_ratio(short, long, name=''):
def get_i(name):
    name = name.split('/')[-1].split('.')[0]
    for a in name.split('_'):
        if a.isdigit():
            return int(a)

def single_bin_mult(short, long, name=None, div=0):
    ''' do [short] (function) [long]
    where function is:
    div == 0: multiply
    div == 1: divide
    div == -1: divide long by short'''
    # assume that the histograms have the same binning and that
    # h_truth contains the range of h_corr

    if name is None:
        name = f'{short.GetName()}x{long.GetName()}_div{div}'

    h_mult = short.Clone(name)
    h_mult.Reset()

    x0 = short.GetXaxis().GetBinLowEdge(1)
    x1 = short.GetXaxis().GetBinUpEdge(short.GetNbinsX())

    for i_long in range(1, long.GetNbinsX()+1):
        x = long.GetXaxis().GetBinCenter(i_long)
        if x < x0 or x > x1:
            continue
        i_short = short.GetXaxis().FindBin(x)
        
        if x != short.GetXaxis().GetBinCenter(i_short):
            print('b0 Error in binning in single_bin_mult')
            sys.exit()

        v_short = short.GetBinContent(i_short)
        v_long  = long.GetBinContent(i_long)
        if v_short == 0 or v_long == 0:
            continue
        e_short = short.GetBinError(i_short)
        e_long  = long.GetBinError(i_long)
        if div == 0:
            v_mult = v_long * v_short
        elif div == 1:
            v_mult = v_short / v_long
        elif div == -1:
            v_mult = v_long / v_short
        else:
            print('b1 Error in division in single_bin_mult')
            sys.exit()
        e_mult = np.sqrt((e_short/v_short)**2 + (e_long/v_long)**2) * v_mult
        h_mult.SetBinContent(i_short, v_mult)
        h_mult.SetBinError(i_short, e_mult)
    return h_mult

def make_err(href, list_others, fillcolor, fillalpha):
    ''' generate histogram with errors that is mean-sum-errors between hred and list_others'''
    herr = href.Clone(f"err_{href.GetName()}")
    sumsq = [0 for x in range(href.GetNbinsX())]
    for hg in list_others:
        for i in range(1, href.GetNbinsX()+1):
            sumsq[i-1] += (href.GetBinContent(i)-hg.GetBinContent(i))**2.
    for i in range(1, href.GetNbinsX()+1):
        herr.SetBinError(i, (sumsq[i-1]/len(list_others))**0.5)
        # herr.SetBinError(i, 0.0001)
    herr.SetMarkerStyle(rt.kDot)
    herr.SetMarkerColorAlpha(rt.kWhite, 1.)
    herr.SetFillColorAlpha(fillcolor, fillalpha)
    return herr

def df_to_th1d_cnt(df, column, bins, name, errs = False, debug=None):
    '''Build a root TH1D form df['column'] in with bins==(nbins, lobin, hibin)'''

    lbins = np.linspace(bins[1], bins[2], bins[0]+1)
    H, _ = np.histogram(df[column], bins=lbins)

    h1 = rt.TH1D(name, f'h1;{column}', *bins)
    h1.Sumw2()

    for (i,), val in np.ndenumerate(H):
        if val == 0:
            continue
        if lbins[i] != h1.GetXaxis().GetBinLowEdge(i+1):
            print('a0 Error in binning')
            sys.exit()
        h1.SetBinContent(i+1, val)
        if debug and (i+1) == debug[2]:
            with open(debug[0], 'a') as fout:
                fout.write(f'{debug[0]} ifile[{debug[1]}] bin {i+1}: {val}\n')
            print(f'ifile[{debug[0]}] bin {i+1}: {val}\n')
    if errs:
        for i in range(1, h1.GetNbinsX()+1):
            h1.SetBinError(i, np.sqrt(h1.GetBinContent(i)))
    return h1

def df_to_th2d_cnt(df, column_x, column_y, bins_x, bins_y, name):
    '''Build a root TH2D form df['column_x']:df['column_y'] in with bins==(nbins, lobin, hibin)'''
    # Get the weighting parameters, too
    lbins_x = np.linspace(bins_x[1], bins_x[2], bins_x[0]+1)
    lbins_y = np.linspace(bins_y[1], bins_y[2], bins_y[0]+1)
    H, _, _ = np.histogram2d(df[column_x], df[column_y], bins=[lbins_x,lbins_y])

    h2 = rt.TH2D(name, f'h2;{column_x};{column_y}', *bins_x, *bins_y)
    h2.Sumw2()

    for (i, j), val in np.ndenumerate(H):
        if val == 0:
            continue
        if (lbins_x[i] != h2.GetXaxis().GetBinLowEdge(i+1)) or (lbins_y[j] != h2.GetYaxis().GetBinLowEdge(j+1)):
            print('a2 Error in binning')
            sys.exit()
        h2.SetBinContent(i+1, j+1, val)

    return h2

def remove_bin(hx, hy, h2, ix, iy):
    '''remove bin ix, iy from h2, and also subtract it from hx and hy'''
    n_remove = h2.GetBinContent(ix, iy)
    h2.SetBinContent(ix, iy, 0)
    hx.SetBinContent(ix, hx.GetBinContent(ix) - n_remove)
    hy.SetBinContent(iy, hy.GetBinContent(iy) - n_remove)

def remove_outliers(df, columns, sigma=4):
    print(f'FIXME df: {df.columns}')

    mean0 = df[columns[0]].mean()
    std0 = df[columns[0]].std()
    lower_bound0 = mean0 - sigma * std0
    upper_bound0 = mean0+ sigma * std0

    cut = (df[columns[0]] >= lower_bound0) & (df[columns[0]] <  upper_bound0)

    for col in columns[1:]:
        mean = df[col].mean()
        std = df[col].std()
        lower_bound = mean - sigma * std
        upper_bound = mean + sigma * std
        cut = cut & (df[col] >= lower_bound) & (df[col] <  upper_bound)

    # mean1 = df[columns[1]].mean()
    # std1 = df[columns[1]].std()
    # lower_bound1 = mean1 - sigma * std1
    # upper_bound1 = mean1+ sigma * std1

    filtered_df = df[cut]
        #  (df[columns[0]] >= lower_bound0) 
    #    & (df[columns[0]] <  upper_bound0)
    #    & (df[columns[1]] >= lower_bound1) 
    #    & (df[columns[1]] <  upper_bound1)]
    print(f'Before {df.shape[0]}, After {filtered_df.shape[0]}, Ratio {filtered_df.shape[0] / df.shape[0]}')
    return filtered_df

def remove_low_stat_bins(hx=None, hy=None, h2=None, min_entry=1, debug=None):
    '''
    Takes in a TH2D histogram
    - Find the bin with the maximum content
    - Working from that row down, for each row:
       - finds the column entry with maximum input
       - sees if that column is within the range of the precvious row
         (if not, delete all following)
       - goes to the left and right until it finds an entry < min_entry
         for then deletes all entries going outward
    - Does the same thing from the maximum cell and working upwards
    '''
    max_content = -1
    max_bin_x = -1
    max_bin_y = -1

    if debug:
        print(f'FIXME:{debug}')
        with open(debug[0], 'a') as fout:
            fout.write(f'file[{debug[1]:2s}] bin[{debug[2]:d}]  pre-clean truth:  {hx.GetBinContent(debug[2]):9.1f} corr: {hy.GetBinContent(debug[2]):9.1f} <\n')

    # Loop through all bins
    for ix in range(1, h2.GetNbinsX() + 1):
        for iy in range(1, h2.GetNbinsY() + 1):
            content = h2.GetBinContent(ix, iy)
            if content > max_content:
                max_content = content
                max_bin_x = ix
                max_bin_y = iy

    if max_content < min_entry:
        for irow in range(1, h2.GetNbinsY()+1):
            for icol in range(1, h2.GetNbinsX()+1):
                remove_bin(hx, hy, h2, icol, irow)
        return 

    # loop down from max_bin_y
    col_first = max_bin_x
    col_last  = max_bin_x

    cut_all = False
    for irow in range(max_bin_y, 0, -1):
        # print(f'row {irow} {cut_all}')
        if cut_all:
            for icol in range(1, h2.GetNbinsX()+1):
                remove_bin(hx, hy, h2, icol, irow)
            continue

        imax = 0
        nmax = 0
        for icol in range(1, h2.GetNbinsY()+1):
            val = h2.GetBinContent(icol, irow)
            if val > min_entry and val > nmax:
                nmax = h2.GetBinContent(icol, irow)
                imax = icol

        if imax == 0 or (not (imax >= col_first) and (imax <= col_last)):
            cut_all = True
            for icol in range(1, h2.GetNbinsX()+1):
                remove_bin(hx, hy, h2, icol, irow)
            continue
        
        # trim ends of rows
        # print(f'row {irow} {imax} nmax: {nmax}   {col_first} {col_last}')
        delete = False
        col_first = imax
        col_last  = imax
        for icol in range(imax, h2.GetNbinsX()+1):
            if delete and h2.GetBinContent(icol, irow) > 0:
                remove_bin(hx, hy, h2, icol, irow)
            elif h2.GetBinContent(icol, irow) < min_entry:
                delete = True
                remove_bin(hx, hy, h2, icol, irow)
            else:
                col_last = icol

        delete = False
        for icol in range(imax, 0, -1):
            if delete and h2.GetBinContent(icol, irow) > 0:
                remove_bin(hx, hy, h2, icol, irow)
            elif h2.GetBinContent(icol, irow) < min_entry:
                delete = True
                remove_bin(hx, hy, h2, icol, irow)
            else:
                col_first = icol

    # loop up from max_bin_y
    col_first = max_bin_x
    col_last  = max_bin_x

    cut_all = False
    for irow in range(max_bin_y, h2.GetNbinsY()+1):
        if cut_all:
            for icol in range(1, h2.GetNbinsX()+1):
                remove_bin(hx, hy, h2, icol, irow)
            continue

        imax = 0
        nmax = 0
        for icol in range(1, h2.GetNbinsY()+1):
            val = h2.GetBinContent(icol, irow)
            if val > min_entry and val > nmax:
                nmax = h2.GetBinContent(icol, irow)
                imax = icol

        if imax == 0 or (not (imax >= col_first) and (imax <= col_last)):
            cut_all = True
            for icol in range(1, h2.GetNbinsX()+1):
                remove_bin(hx, hy, h2, icol, irow)
            continue
        
        # trim ends of rows
        # print(f'row {irow} {imax} nmax: {nmax}   {col_first} {col_last}')
        delete = False
        col_first = imax
        col_last  = imax
        for icol in range(imax, h2.GetNbinsX()+1):
            if delete and h2.GetBinContent(icol, irow) > 0:
                remove_bin(hx, hy, h2, icol, irow)
            elif h2.GetBinContent(icol, irow) < min_entry:
                delete = True
                remove_bin(hx, hy, h2, icol, irow)
            else:
                col_last = icol

        delete = False
        for icol in range(imax, 0, -1):
            if delete and h2.GetBinContent(icol, irow) > 0:
                remove_bin(hx, hy, h2, icol, irow)
            elif h2.GetBinContent(icol, irow) < min_entry:
                delete = True
                remove_bin(hx, hy, h2, icol, irow)
            else:
                col_first = icol
    if debug:
        with open(debug[0], 'a') as fout:
            fout.write(f'file[{debug[1]:2s}] bin[{debug[2]:d}]  post-clean truth: {hx.GetBinContent(debug[2]):9.1f} corr: {hy.GetBinContent(debug[2]):9.1f} >\n')

def style_hg(hg, marker_style, color, write=False, alpha=None):
    '''Set marker-style and marker-color of TH1D; 
    also set range to 10% below lowest and 10% above highest bin content (excepting 0)'''
    hg.SetMarkerStyle(marker_style)

    if alpha:
        hg.SetMarkerColorAlpha(color,alpha)
        hg.SetLineColorAlpha(color,alpha)
    else:
        hg.SetMarkerColor(color)
        hg.SetLineColor(color)

    max_val = hg.GetMaximum()
    min_val = max_val

    if min_val == 0: 
        return
    for i in range(1, hg.GetNbinsX() + 1):
        val = hg.GetBinContent(i)
        if val == 0:
            continue
        if val < min_val:
            min_val = val

    hg.GetYaxis().SetRangeUser(0.1*min_val, 1.1*max_val)
    if write:
        hg.Write()

def err_hist(h_corr, h_truth, h_resp):
    '''Set the errors in TH1D*, TH1D*, TH2D* input'''
    for irow in range(1, h_resp.GetNbinsY()+1):
        c = h_corr.GetBinContent(irow)
        if c == 0:
            continue
        if c < 0:
            print(('h_corr neg: ',c))
            continue
        h_corr.SetBinError(irow, np.sqrt(c))
    for icol in range(1, h_resp.GetNbinsX()+1):
        c = h_truth.GetBinContent(icol)
        if c == 0:
            continue
        h_truth.SetBinError(icol, np.sqrt(c))
    for irow in range(1, h_resp.GetNbinsY()+1):
        for icol in range(1, h_resp.GetNbinsX()+1):
            c = h_resp.GetBinContent(icol, irow)
            if c == 0:
                continue
            h_resp.SetBinError(icol, irow, np.sqrt(c))

def weight_hist(h_corr, h_truth, h_resp, weight, debug=False):
    '''multipy the weights to TH1*, TH1*, TH2* the scale_file json 'scale_factor' '''
    if debug:
        print(f"weight factor: {weight:.6g}")
        try:
            h_corr_int = h_corr.Integral()
            h_truth_int = h_truth.Integral()
            h_resp_int = h_resp.Integral()
            # print(f'\nints before: {h_corr_int:.4g} {h_truth_int:.4g} {h_resp_int:.4g}')
            # print(f'scaling h_corr before: {h_corr.GetName()} {h_corr.Integral()}')
        except:
            pass
    h_corr.Scale(weight)

    # print(f'scaling h_corr after: {h_corr.GetName()} {h_corr.Integral()}')
    h_truth.Scale(weight)
    h_resp.Scale(weight)
    if debug:
        try:
            print(f'ratio ints after: {h_corr.Integral()/h_corr_int} {h_truth.Integral()/h_truth_int} {h_resp.Integral()/h_resp_int}')
        except:
            pass
    return weight


def BuildResponse(bins_x=[120,0,120], bins_y=[120,0,120], Bratio=0.4, min_bin_content=10, truth_pt_range=(0,1000.), corr_pt_range=(0,1000.),  true_corr_dir='.', root_ofile='response.root',save_all=False, outlier_sigma=4., draw_stacked=True, comp_cross_file=None, min_rhoApt=0., clean="_clean"):
    # main program
    '''
    Fill a 1D histogram of Truth (no bounds on Corr)
    Fill a 1D histogram of Corr  (no bounds on Truth)
    Fill a 2D histogram of Corr vs Truth

    Remove bins in 2D will less than bin content and then subtract 
    from the 1D histograms
    '''

    if comp_cross_file:
        with open(comp_cross_file, 'r') as f:
            comp_cross = json.load(f)
    else:
        comp_cross = None

    print(f'FIXME comp_cross: {comp_cross}')


    fout = rt.TFile(root_ofile, 'recreate')
    print(f'Building {root_ofile}')

    if not os.path.isdir(true_corr_dir):
        print(f'Cannot find directory {true_corr_dir}')
        sys.exit()

    if draw_stacked:
        df_cumulative = pd.DataFrame()

    first = True
    print(f'{true_corr_dir}/truth_reco_corr/data_*{clean}.parquet')

    sum_cross_section = 0

    for file in sorted(
        glob(f'{true_corr_dir}/truth_reco_corr/data_*{clean}.parquet'), key=lambda x: get_i(x)):

        ifile = get_i(file)
        print(f'This file: {file} and i: {ifile}')
        if int(ifile) > 14:
            continue


        df = pd.read_parquet(file)
        print(f'FIXME A0: number now {df.shape[0]}')
        # print(f'df columns: {df.columns}')
        df = remove_outliers(df, ('truth_pt', 'reco_pt'), outlier_sigma)
        print(f'FIXME A1: number now {df.shape[0]}')

        cuts = df['truth_pt'] > -100. # this is hack to get the cuts always true
        if truth_pt_range:
            cuts = cuts & (df['truth_pt'] > truth_pt_range[0]) & (df['truth_pt'] < truth_pt_range[1])
        if corr_pt_range:
            cuts = cuts & (df['corr_pt'] > corr_pt_range[0]) & (df['corr_pt'] < corr_pt_range[1])
        if min_rhoApt:
            cuts = cuts & ( (df['reco_pt']-df['rho_bkg_thermal']*df['reco_area']) > min_rhoApt)
        
        df = df[cuts].copy()
        print(f'FIXME A2: number now {df.shape[0]}')

        h_corr  = df_to_th1d_cnt(df, 'corr_pt', bins_x, f'corr_pt_{ifile}', debug=("h_corr_bin_9", ifile, 9))
        h_truth = df_to_th1d_cnt(df, 'truth_pt', bins_y, f'truth_pt_{ifile}', debug=("h_truth_bin_9",  ifile, 9))
        h_resp= df_to_th2d_cnt(df, 'corr_pt', 'truth_pt', bins_x, bins_y, f'response_{ifile}')

        A_df = df.sample(frac=Bratio)
        A_corr  = df_to_th1d_cnt(A_df, 'corr_pt', bins_x, f'A_corr_pt_{ifile}')
        A_truth = df_to_th1d_cnt(A_df, 'truth_pt', bins_y, f'A_truthpt_{ifile}')
        A_resp= df_to_th2d_cnt(A_df, 'corr_pt', 'truth_pt', bins_x, bins_y, f'A_response_{ifile}')

        B_df = df.drop(A_df.index)
        B_corr  = df_to_th1d_cnt(B_df, 'corr_pt', bins_x, f'B_corr_pt_{ifile}')
        B_truth = df_to_th1d_cnt(B_df, 'truth_pt', bins_y, f'B_truthpt_{ifile}')
        B_resp= df_to_th2d_cnt(B_df, 'corr_pt', 'truth_pt', bins_x, bins_y, f'B_response_{ifile}')

        # clean: remove bins with less than min_bin_content in h_resp and
        # propagate that loss to the h_truth and h_corr
        if (min_bin_content > 0):
            remove_low_stat_bins(h_corr, h_truth, h_resp, min_bin_content, debug=("h_rm_bin_9", ifile, 9))
            remove_low_stat_bins(A_corr, A_truth, A_resp, min_bin_content)
            remove_low_stat_bins(B_corr, B_truth, B_resp, min_bin_content)

        # add sqrt to the histograms
        err_hist(h_corr, h_truth, h_resp)
        err_hist(A_corr, A_truth, A_resp)
        err_hist(B_corr, B_truth, B_resp)

        # scale the histograms
        with open(f'{true_corr_dir}/../Xsec_groups/scale_data_{ifile}.json', 'r') as f:
            weight = json.load(f)['scale_factor']
        if comp_cross:
            weight *= comp_cross[f'file_{ifile}']
        weight_hist(h_corr, h_truth, h_resp, weight, debug=False)
        weight_hist(A_corr, A_truth, A_resp, weight, debug=False)
        weight_hist(B_corr, B_truth, B_resp, weight, debug=False)

        sum_cross_section += weight * h_truth.Integral()

        if draw_stacked:
            print()
            df['weight'] = w
            print(f'weight {w} file {file}')
            df_cumulative = pd.concat([df_cumulative, df])
            

        if first:
            all_h_truth = h_truth.Clone('truth_pt')
            all_h_corr  = h_corr.Clone('corr_pt')
            all_h_resp = h_resp.Clone('response')

            all_A_truth = A_truth.Clone('A_truth_pt')
            all_A_corr  = A_corr.Clone('A_corr_pt')
            all_A_resp = A_resp.Clone('A_response')

            all_B_truth = B_truth.Clone('B_truth_pt')
            all_B_corr  = B_corr.Clone('B_corr_pt')
            all_B_resp = B_resp.Clone('B_response')
            first = False
        else:
            all_h_truth.Add(h_truth)
            all_h_corr.Add(h_corr)
            all_h_resp.Add(h_resp)

            all_A_truth.Add(A_truth)
            all_A_corr.Add(A_corr)
            all_A_resp.Add(A_resp)

            all_B_truth.Add(B_truth)
            all_B_corr.Add(B_corr)
            all_B_resp.Add(B_resp)

        if save_all:
            h_truth.Write()
            h_corr.Write()
            h_resp.Write()

            A_truth.Write()
            A_corr.Write()
            A_resp.Write()

            B_truth.Write()
            B_corr.Write()
            B_resp.Write()

    if draw_stacked:
        import matplotlib.pyplot as plt
        bins = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50]
        labels = [f'{bins[i]}-{bins[i+1]}' for i in range(len(bins)-1)]
        df = df_cumulative
        # Create a new column to categorize truth_pt into bins
        df['truth_pt_bin'] = pd.cut(df['truth_pt'], bins=bins, labels=labels, right=False)

        # Group the DataFrame by the new column
        grouped = df.groupby('truth_pt_bin')

        # Prepare data for the stacked histogram
        hist_data = []
        weights = []
        print(f'FIXME {np.max(df["weight"])} and {np.min(df["weight"])}')
        for name, group in grouped:
            print(f'group {name} {group.shape}')
            print('group pt')
            print(group['corr_pt'])
            print('weight')
            print(np.max(group['weight']))
            print(np.min(group['weight']))
            print(group['weight'])
            hist_data.append(group['corr_pt'])
            weights.append(group['weight'])

        # Plot the stacked histogram
        plt.figure(figsize=(10, 6))
        plt.hist(hist_data, bins=50, stacked=True, weights=weights, label=labels, edgecolor='black')

        # Add labels and title
        # plt.yscale('log')
        plt.xlabel('pt_rhoA')
        plt.ylabel('Cross Section')
        plt.xlabel(r'ML corrected pT')
        plt.title('Stacked Histogram of pTcorr Weighted by weight')
        plt.legend(title='truth_pt bins')

        # Show the plot
        plt.yscale('log')
        plt.savefig(f'{root_ofile.replace(".root", "_stacked.png")}')
        plt.show()


    style_hg(all_h_truth, rt.kOpenCircle, rt.kBlack)
    style_hg(all_A_truth, rt.kOpenCircle, rt.kRed)
    style_hg(all_B_truth, rt.kOpenCircle, rt.kBlue)

    style_hg(all_h_corr, rt.kOpenSquare, rt.kBlack)
    style_hg(all_A_corr, rt.kOpenSquare, rt.kRed)
    style_hg(all_B_corr, rt.kOpenSquare, rt.kBlue)

    print(f' PEAR FIXME truth: {all_h_truth.GetBinError(3)} corr: {all_h_corr.GetBinError(3)} resp: {all_h_resp.GetBinError(3)}')
    all_h_truth.Write()
    all_h_corr.Write()
    all_h_resp.Write()

    all_A_truth.Write()
    all_A_corr.Write()
    all_A_resp.Write()

    all_B_truth.Write()
    all_B_corr.Write()
    all_B_resp.Write()

    fout.Save()
    print(f'Done writing {fout}')

def create_RAA_like(pp_file, AuAu_file, fout_name, prefix='', niter=4, proj_truth=False):
    ''' Write out:
    truth_pp
    sbe_pp
    ruu_response_pp

    corr_pp_sbe
    corr_pp_ruu2
    corr_pp_ruu4
    corr_pp_ruu6

    rat_corr_pp_sbe
    rat_corr_pp_ruu2
    rat_corr_pp_ruu4
    rat_corr_pp_ruu6

    truth_AuAu
    corr_AuAu_sbe
    corr_AuAu_ruu2
    corr_AuAu_ruu4
    corr_AuAu_ruu6

    rat_truth_AuAu
    rat_corr_AuAu_sbe
    rat_corr_AuAu_ruu2
    rat_corr_AuAu_ruu4
    rat_corr_AuAu_ruu6
    '''

    f_pp   = rt.TFile.Open(pp_file)
    _T_pp = f_pp.Get(f"{prefix}truth_pt")
    _C_pp = f_pp.Get(f"{prefix}corr_pt")
    _R_pp = f_pp.Get(f"{prefix}response")


    fout = rt.TFile(fout_name, 'recreate')
    fout.cd()
    if proj_truth:
        T_pp = _R_pp.ProjectionY('truth_pp')
    else:
        T_pp = _T_pp.Clone('truth_pp')

    T_pp.Write()
    C_pp = _C_pp.Clone('corr_pp')
    C_pp.Write()
    print(f'FIXME: {_R_pp.GetName()}')
    R_pp = _R_pp.Clone('response_pp')
    print(f'FIXME: {R_pp.GetName()} and Int: {R_pp.Integral()}')
    R_pp.Write()


    f_AuAu = rt.TFile.Open(AuAu_file)
    _T_AuAu = f_AuAu.Get(f"{prefix}truth_pt")
    _C_AuAu = f_AuAu.Get(f"{prefix}corr_pt")
    _R_AuAu = f_AuAu.Get(f"{prefix}response")

    fout.cd()
    T_AuAu = _T_AuAu.Clone('truth_AuAu')
    if proj_truth:
        T_AuAu = _R_AuAu.ProjectionY('truth_AuAu')
    T_AuAu.Write()
    C_AuAu = _C_AuAu.Clone('corr_AuAu')
    C_AuAu.Write()
    print(f'FIXME: {_R_AuAu.GetName()}')
    R_AuAu = _R_AuAu.Clone('response_AuAu')
    print(f'FIXME: {R_AuAu.GetName()} and Int: {R_AuAu.Integral()}')
    R_AuAu.Write()

    sbe_pp = single_bin_mult(C_pp, T_pp, 'sbe_pp', -1)
    style_hg(sbe_pp, rt.kOpenSquare, rt.kBlack, write=True)

    ruu_response_pp_2 = rt.RooUnfoldResponse(C_pp, T_pp, R_pp)
    ruu_response_pp_2.SetName('ruu_response_pp_4')

    ruu_response_pp_4 = rt.RooUnfoldResponse(C_pp, T_pp, R_pp)
    ruu_response_pp_4.SetName('ruu_response_pp_6')

    ruu_response_pp_6 = rt.RooUnfoldResponse(C_pp, T_pp, R_pp)
    ruu_response_pp_6.SetName('ruu_response_pp_8')
    # ruu_response_pp.Write()

    corr_pp_sbe = single_bin_mult(C_pp, sbe_pp, 'corr_pp_sbe')
    style_hg(corr_pp_sbe, rt.kOpenTriangleUp, rt.kBlue, write=True)

    bayes_pp_ruu2 = rt.RooUnfoldBayes(ruu_response_pp_2, C_pp, niter-2)
    corr_pp_ruu2 = bayes_pp_ruu2.Hunfold()
    corr_pp_ruu2.SetName("corr_pp_ruu_2")
    style_hg(corr_pp_ruu2, rt.kOpenTriangleDown, rt.kRed, write=True)

    bayes_pp_ruu4 = rt.RooUnfoldBayes(ruu_response_pp_4, C_pp, niter)
    corr_pp_ruu4 = bayes_pp_ruu4.Hunfold()
    corr_pp_ruu4.SetName("corr_pp_ruu_4")
    style_hg(corr_pp_ruu4, rt.kOpenStar, (rt.kGreen+2), write=True)

    bayes_pp_ruu6 = rt.RooUnfoldBayes(ruu_response_pp_6, C_pp, niter+2)
    corr_pp_ruu6 = bayes_pp_ruu6.Hunfold()
    corr_pp_ruu6.SetName("corr_pp_ruu_6")
    style_hg(corr_pp_ruu6, rt.kOpenStar, (rt.kMagenta+2), write=True)

    ruu_CHECK = rt.RooUnfoldResponse(C_AuAu, T_AuAu, R_AuAu)
    ruu_CHECK.SetName('ruu_CHECK')
    bayes_CHECK = rt.RooUnfoldBayes(ruu_CHECK, C_AuAu, niter)
    corr_CHECK_ruu4 = bayes_CHECK.Hunfold()
    corr_CHECK_ruu4.SetName("corr_AuCHECK_ruu_4")
    style_hg(corr_CHECK_ruu4, rt.kOpenStar, (rt.kGreen+2), write=True)


    rat_corr_pp_sbe = single_bin_mult(corr_pp_sbe, T_pp, 'rat_corr_pp_sbe', 1)
    rat_corr_pp_sbe.Write()

    rat_corr_pp_ruu2 = single_bin_mult(corr_pp_ruu2, T_pp, 'rat_corr_pp_ruu2', div=1)
    rat_corr_pp_ruu2.GetYaxis().SetRangeUser(0., 1.2)
    rat_corr_pp_ruu2.Write()

    rat_corr_pp_ruu4 = single_bin_mult(corr_pp_ruu4, T_pp, 'rat_corr_pp_ruu4', div=1)
    rat_corr_pp_ruu4.GetYaxis().SetRangeUser(0., 1.2)
    rat_corr_pp_ruu4.Write()

    rat_corr_pp_ruu6 = single_bin_mult(corr_pp_ruu6, T_pp, 'rat_corr_pp_ruu6', div=1)
    rat_corr_pp_ruu6.GetYaxis().SetRangeUser(0., 1.2)
    rat_corr_pp_ruu6.Write()

    corr_AuAu_sbe = single_bin_mult(C_AuAu, sbe_pp, 'corr_AuAu_sbe')
    print(corr_AuAu_sbe.GetName())
    style_hg(corr_AuAu_sbe, rt.kOpenTriangleUp, rt.kBlue, write=True)
    corr_AuAu_sbe.Write()

    bayes_AuAu_ruu2 = rt.RooUnfoldBayes(ruu_response_pp_2, C_AuAu, niter-2)
    corr_AuAu_ruu2 = bayes_AuAu_ruu2.Hunfold()
    corr_AuAu_ruu2.SetName("corr_AuAu_ruu_2")
    style_hg(corr_AuAu_ruu2, rt.kOpenCircle, rt.kRed, write=True)
    corr_AuAu_ruu2.Write()

    bayes_AuAu_ruu4 = rt.RooUnfoldBayes(ruu_response_pp_4,C_AuAu,  niter)
    corr_AuAu_ruu4 = bayes_AuAu_ruu4.Hunfold()
    corr_AuAu_ruu4.SetName("corr_AuAu_ruu_4")
    style_hg(corr_AuAu_ruu4, rt.kOpenCircle, rt.kGreen+2, write=True)
    corr_AuAu_ruu4.Write()

    bayes_AuAu_ruu6 = rt.RooUnfoldBayes(ruu_response_pp_6, C_AuAu, niter+2)
    corr_AuAu_ruu6 = bayes_AuAu_ruu6.Hunfold()
    corr_AuAu_ruu6.SetName("corr_AuAu_ruu_6")
    style_hg(corr_AuAu_ruu6, rt.kOpenCircle, rt.kOrange+2, write=True)
    corr_AuAu_ruu6.Write()

    bayes_AuAu_error = make_err(corr_AuAu_ruu4, [corr_AuAu_ruu2, corr_AuAu_ruu6], rt.kGreen+2, 0.5)
    bayes_AuAu_error.SetName('bayes_AuAu_error')
    bayes_AuAu_error.Write()

    rat_truth_AuAu = single_bin_mult(_T_AuAu, _T_pp, 'rat_truth_AuAu', div=1)
    style_hg(rat_truth_AuAu, rt.kOpenSquare, rt.kRed, write=True) 
    rat_truth_AuAu.GetYaxis().SetRangeUser(0., 1.2)
    rat_truth_AuAu.Write()

    rat_corr_AuAu_sbe = single_bin_mult(corr_AuAu_sbe, T_pp, 'rat_corr_AuAu_sbe',div=1)
    style_hg(rat_corr_AuAu_sbe, rt.kOpenTriangleUp, rt.kBlue, write=True)
    rat_corr_AuAu_sbe.GetYaxis().SetRangeUser(0., 1.2)
    rat_corr_AuAu_sbe.Write()

    rat_corr_AuAu_ruu2 = single_bin_mult(corr_AuAu_ruu2, T_pp, 'rat_corr_AuAu_ruu2',div=1)
    style_hg(rat_corr_AuAu_ruu2, rt.kOpenCircle, rt.kGreen, write=True)
    rat_corr_AuAu_ruu2.GetYaxis().SetRangeUser(0., 1.2)
    rat_corr_AuAu_ruu2.Write()

    rat_corr_AuAu_ruu4 = single_bin_mult(corr_AuAu_ruu4, T_pp, 'rat_corr_AuAu_ruu4',div=1)
    style_hg(rat_corr_AuAu_ruu4, rt.kOpenCircle, rt.kGreen+2, write=True)
    rat_corr_AuAu_ruu4.GetYaxis().SetRangeUser(0., 1.2)
    rat_corr_AuAu_ruu4.Write()

    rat_corr_AuAu_ruu6 = single_bin_mult(corr_AuAu_ruu6, T_pp, 'rat_corr_AuAu_ruu6',div=1)
    style_hg(rat_corr_AuAu_ruu6, rt.kOpenCircle, rt.kOrange+2, write=True)
    rat_corr_AuAu_ruu6.GetYaxis().SetRangeUser(0., 1.2)
    rat_corr_AuAu_ruu6.Write()

    rat_corr_AuAu_ruuCHECK = single_bin_mult(corr_CHECK_ruu4, T_pp, 'rat_CHECK_AuAu_ruu4',div=1)
    style_hg(rat_corr_AuAu_ruuCHECK, rt.kOpenCircle, rt.kOrange+2, write=True)
    rat_corr_AuAu_ruuCHECK.GetYaxis().SetRangeUser(0., 1.2)
    rat_corr_AuAu_ruuCHECK.Write()

    err_rat_bayes_AuAu = make_err(rat_corr_AuAu_ruu4, [rat_corr_AuAu_ruu2, rat_corr_AuAu_ruu6], rt.kGreen+2, 0.5)
    err_rat_bayes_AuAu.SetName('err_rat_bayes_AuAu')
    err_rat_bayes_AuAu.GetYaxis().SetRangeUser(0., 1.2)
    err_rat_bayes_AuAu.Write()

    fout.Save()

def draw_th1d(hg, marker, color, size, legend=None, text='', alpha=1., first=False, title=''):
    hg.SetStats(0)
    hg.SetMarkerStyle(marker)
    hg.SetMarkerColorAlpha(color,alpha)
    hg.SetLineColorAlpha(color,alpha)
    hg.SetMarkerSize(size)
    hg.Draw('PE' if first else 'PE SAME')
    hg.SetTitle(title)
    if legend:
        legend.AddEntry(hg, text, 'p')

def scale_semilogy(hg, low=None, high=None, debug=None):
    vmin = 0
    vmax = 0
    for i in range(1, hg.GetNbinsX()+1):
        val = hg.GetBinContent(i)
        if val == 0:
            continue
        if vmin == 0:
            vmin = val
            vmax = val
            continue
        if val < vmin:
            vmin = val
        if val > vmax:
            vmax = val
    if not low:
        low = vmin*0.5
    if not high:
        high = vmax*1.1
    hg.GetYaxis().SetRangeUser(low,high)
    if debug:
        print(f'low: {low} high: {high}')
    return low, high

def scale_yaxis(hg, low=None, high=None, debug=None):
    vmin = hg.GetBinContent(1)
    vmax = hg.GetBinContent(1)
    for i in range(1, hg.GetNbinsX()+1):
        val = hg.GetBinContent(i)
        if val < vmin:
            vmin = val
        if val > vmax:
            vmax = val
    if not low:
        low = vmin*0.5
    if not high:
        high = vmax*1.1
    hg.GetYaxis().SetRangeUser(low,high)
    if debug:
        print(f'low: {low} high: {high}')
    return low, high

def get_stats(which, lowpt=45, highpt=50, max_ifile=20, clean="_clean"):
    if which=='pp':
        dir_name = '../NoBrickInput/M_AB_method_NoBrick'
    else:
        dir_name = '../Threep5FermiInput/M_AB_method_NoBrick'
    results = dict()
    for file in glob(f'{dir_name}/truth_reco_corr/data_*{clean}.parquet'):
        ifile = get_i(file)
        if ifile > max_ifile:
            continue
        df = pd.read_parquet(file)
        cnt_truth = df[(df['truth_pt']>=lowpt) & (df['truth_pt']<highpt)].shape[0]
        cnt_corr  = df[(df['corr_pt']>=lowpt) & (df['corr_pt']<highpt)].shape[0]
        with open(f'{dir_name}/../Xsec_groups/scale_data_{ifile}.json', 'r') as f:
            j = json.load(f)
            weights = j['scale_factor']
        if (cnt_truth != 0) and (cnt_corr != 0):
            rel_err = (1./cnt_truth + 1./cnt_corr)**0.5
            sbe = cnt_truth/cnt_corr
        else:
            rel_err = 0.
            sbe = 0.
        if cnt_truth == 0:
            rerr_truth = 0.
        else:
            rerr_truth = 1./np.sqrt(cnt_truth)
        if cnt_corr == 0:
            rerr_corr = 0.
        else:
            rerr_corr = 1./np.sqrt(cnt_corr)
        results[ifile]={'W':weights, 'cnt_truth':cnt_truth,'rerr_truth':rerr_truth, 'cnt_corr':cnt_corr, 'rerr_corr':rerr_corr, 'SBE':sbe, 'rel_err':rel_err}

    cnt_truth = 0
    e2_truth = 0
    cnt_corr = 0
    e2_corr = 0
    for ifile in sorted(results.keys()):
        W = results[ifile]['W']
        T = results[ifile]['cnt_truth']
        eT = np.sqrt(T)
        C = results[ifile]['cnt_corr']
        eC = np.sqrt(C)

        cnt_truth += W*T
        e2_truth += W**2*eT**2
        if False:
            print(f'{ifile} e2_truth: {e2_truth}  rel_err: {np.sqrt(e2_truth)/cnt_truth} W: {W} T: {cnt_truth} eT: {eT}')
        cnt_corr += W*C
        e2_corr += W**2*eC**2

    rerr_truth = np.sqrt(e2_truth)/cnt_truth
    rerr_corr = np.sqrt(e2_corr)/cnt_corr
    results[100] = {'cnt_truth':cnt_truth, 'rerr_truth':rerr_truth, 'cnt_corr':cnt_corr, 'rerr_corr':rerr_corr, 'SBE':cnt_truth/cnt_corr, 'rel_err':(rerr_truth**2+rerr_corr**2)**0.5, 'W':0.}
    return results

def print_stats(stats, ilow=0, ihi=60, tag=''):
    print(f'{"":14s} {"file":5s} {"Xsec-fac":10s} {"n_truth":30s} {"n_corr":30s} {"SBE":10s} {"rel_err":10s}')
    for ifile in range(ilow, ihi+1):
        if ifile not in stats:
            continue
        S=stats[ifile]
        # print stats[ifile]
        print(f'{tag:14s} {ifile:5d} {S["W"]:10.4g} {S["cnt_truth"]:10d}[w: {(S["cnt_truth"]*S["W"]):13.8g}] {S["cnt_corr"]:10d}[w: {(S["cnt_corr"]*S["W"]):13.8g}] {S["SBE"]:10.5f} {S["rel_err"]:10.4g}')
    ifile = 100
    S = stats[ifile]
    print('HERE')
    print(S)
    print(f'{tag:14s} {ifile:5d} {S["W"]:10.4g} {S["cnt_truth"]:10d}[w: {(S["cnt_truth"]*S["W"]):13.8g}] {S["cnt_corr"]:10d}[w: {(S["cnt_corr"]*S["W"]):13.8g}] {S["SBE"]:10.5f} {S["rel_err"]:10.4g}')
    # print(f'{ifile:5d} {S["W"]:10.4g} {S["cnt_truth"]:10.4g} {S["cnt_corr"]:10.4g} {S["SBE"]:10.5f} {S["rel_err"]:10.4g}')

def plot_SBE_ratios(dir_name, tag, outfile, details=''):
    pp_file   = rt.TFile.Open(f'{dir_name}/pp_{tag}.root')
    AuAu_file = rt.TFile.Open(f'{dir_name}/AuAu_{tag}.root')

    pp_truth   = pp_file.Get('truth_pt')
    pp_corr    = pp_file.Get('corr_pt')

    AuAu_truth = AuAu_file.Get('truth_pt')
    AuAu_corr  = AuAu_file.Get('corr_pt')

    alpha = 0.6
    # shoe-horn in roounfold as well
    pp_R = pp_file.Get('response')
    ruuUnfResp = rt.RooUnfoldResponse(pp_corr, pp_truth, pp_R)
    ruu_5 = rt.RooUnfoldBayes(ruuUnfResp, AuAu_corr, 5)
    hg_5 = ruu_5.Hunfold()
    hg_5.SetName('corr_pp_ruu_5')
    rat_5 = single_bin_mult(hg_5, AuAu_truth, 'rat_corr_pp_ruu_5', div=1)
    style_hg(hg_5,  rt.kFullStar, rt.kGreen+2, 0.8, alpha=alpha)
    style_hg(rat_5, rt.kFullStar, rt.kGreen+2, 0.8, alpha=alpha)

    ruu_15 = rt.RooUnfoldBayes(ruuUnfResp, AuAu_corr, 15)
    hg_15 = ruu_15.Hunfold()
    hg_15.SetName('corr_pp_ruu_15')
    rat_15 = single_bin_mult(hg_15, AuAu_truth, 'rat_corr_pp_ruu_15', div=1)
    style_hg(hg_15,  rt.kFullStar, rt.kGray+2, 0.8, alpha=alpha)
    style_hg(rat_15, rt.kFullStar, rt.kGray+2, 0.8, alpha=alpha)

    ruu_50 = rt.RooUnfoldBayes(ruuUnfResp, AuAu_corr, 50)
    hg_50 = ruu_50.Hunfold()
    hg_50.SetName('corr_pp_ruu_50')
    rat_50 = single_bin_mult(hg_50, AuAu_truth, 'rat_corr_pp_ruu_50', div=1)
    style_hg(hg_50,  rt.kFullStar, rt.kMagenta+2, 0.8, alpha=alpha)
    style_hg(rat_50, rt.kFullStar, rt.kMagenta+2, 0.8, alpha=alpha)

    max_y = 2*np.max([pp_truth.GetMaximum(), pp_corr.GetMaximum(), AuAu_truth.GetMaximum(), AuAu_corr.GetMaximum()])
    min_y = 0.5*np.min([pp_truth.GetMinimum(), pp_corr.GetMinimum(), AuAu_truth.GetMinimum(), AuAu_corr.GetMinimum()])

    c = rt.TCanvas("c", "Canvas", 600, 600)
    c.SetLeftMargin(0.15)
    c.Draw()

    pad1 = rt.TPad("pad1", "Pad 1", 0, 0.5, 1, 1.0)
    pad1.SetBottomMargin(0)  # Upper and lower plot are joined
    pad1.SetLeftMargin(0.15)
    pad1.Draw()

    legend = rt.TLegend(0.6, 0.6, 0.85, 0.85)

    # 1. Divide the TCanvas into 3 TPads
    pad1 = rt.TPad("pad1", "Pad 1", 0, 0.66, 1, 1.0)
    pad1.SetBottomMargin(0)  # Upper and lower plot are joined
    pad1.SetLeftMargin(0.15)
    pad1.Draw()
    pad1.cd()
    pad1.SetLogy()

    pp_truth.GetYaxis().SetRangeUser(min_y, max_y)
    pp_truth.GetYaxis().SetTitle("#frac{dN_{LeadJet}}{d#it{p}_{T}}")
    pp_truth.GetYaxis().SetLabelSize(0.07)
    pp_truth.GetYaxis().SetTitleSize(0.000)
    pp_truth.GetYaxis().SetTitleOffset(0.130)

    draw_th1d(pp_truth, rt.kOpenCircle, rt.kBlack, 1.0, first=True, legend=legend,text=r'$pp$',title=details)

    draw_th1d(pp_corr, rt.kOpenSquare, rt.kBlack, 1.0 , legend=legend, text='pp NN-corrected')

    draw_th1d(AuAu_truth, rt.kOpenCircle, rt.kRed, 1.0, legend=legend, text='Quenched truth')

    draw_th1d(AuAu_corr, rt.kOpenSquare, rt.kRed, 1.0, legend=legend, text='Quenched NN-corrected')

    hg_5.Draw("PE same")
    hg_15.Draw("PE same")
    hg_50.Draw("PE same")

    legend.Draw()

    leg_stars = rt.TLegend(0.17, 0.07, .4, .3)
    leg_stars.AddEntry(hg_5, 'RooUnfold 5', 'p')
    leg_stars.AddEntry(hg_15, 'RooUnfold 15', 'p')
    leg_stars.AddEntry(hg_50, 'RooUnfold 50', 'p')
    leg_stars.Draw()

    # 2. Divide the TCanvas into 3 TPads
    c.cd()
    pad2 = rt.TPad("pad2", "Pad 2", 0, 0.33, 1, 0.66)
    pad2.SetBottomMargin(0)  # Upper and lower plot are joined
    pad2.SetTopMargin(0)
    pad2.SetLeftMargin(0.15)
    pad2.Draw()
    pad2.SetLogy()
    pad2.cd()

    leg2 = rt.TLegend(0.5, 0.6, 0.85, 0.85)

    sbe_pp = single_bin_mult(pp_corr, pp_truth, 'sbe_pp', -1)
    sbe_AuAu = single_bin_mult(AuAu_corr, AuAu_truth, 'sbe_AuAu', -1)

    lo0, hi0 = scale_semilogy(sbe_pp, 0., debug=False)
    lo1, hi1 = scale_semilogy(sbe_AuAu, 0., debug=False)
    sbe_pp.GetYaxis().SetRangeUser(min(lo0, lo1), max(hi0, hi1))

    # print('FIXME A0')
    # for i in range(1, sbe_pp.GetNbinsX()+1):
        # print(f' sbe_pp {i} {sbe_pp.GetXaxis().GetBinCenter(i)}: {sbe_pp.GetBinContent(i)}')
    # max_y = 1.1*np.max([sbe_pp.GetMaximum(), sbe_AuAu.GetMaximum()])
    # min_y = 0.9*np.min([sbe_pp.GetMinimum(), sbe_AuAu.GetMinimum()])

    sbe_pp.GetYaxis().SetTitle("1/SBE: Truth/Corrected")
    sbe_pp.GetYaxis().SetLabelSize(0.07)
    sbe_pp.GetYaxis().SetTitleSize(0.07)
    draw_th1d(sbe_pp, rt.kOpenCircle, rt.kBlack, 1.0, first=True,
              legend=leg2, text='pp: Corr/Truth')
    sbe_pp.Draw("PE")

    sbe_AuAu.GetYaxis().SetTitle("Ratio: Truth/Corrected")
    draw_th1d(sbe_AuAu, rt.kOpenCircle, rt.kRed, 1.0,
             legend=leg2, text='Quenched: Corr/Truth' )
    # leg2.Draw()

    # 3. Ratio of sbe
    c.cd()
    pad3 = rt.TPad("pad3", "Pad 3", 0, 0.0, 1., 0.33)
    pad3.SetTopMargin(0)
    pad3.SetLeftMargin(0.15)
    pad3.SetBottomMargin(0.25)
    pad3.Draw()
    pad3.cd()

    leg3 = rt.TLegend(0.34, 0.6, 0.60, 0.85)

    rat_sbe = single_bin_mult(sbe_pp, sbe_AuAu, 'rat_sbe', div=1)
    scale_yaxis(rat_sbe, 0.)
    # print(f'pp: {rat_sbe.GetMinimum()} {rat_sbe.GetMaximum()}')
    # for i in range(1, rat_sbe.GetNbinsX()+1):
        # print(f'{i} {rat_sbe.GetXaxis().GetBinCenter(i)}: {rat_sbe.GetBinContent(i)}')
    # rat_sbe.GetYaxis().SetRangeUser(0., 1.1)
    rat_sbe.GetYaxis().SetTitle("(1/SBE):(pp/Quench) stars:ruu/truth")
    rat_sbe.GetXaxis().SetTitle("Jet pT")
    rat_sbe.GetYaxis().SetLabelSize(0.07)
    rat_sbe.GetYaxis().SetTitleSize(0.07)
    rat_sbe.GetXaxis().SetLabelSize(0.07)
    rat_sbe.GetXaxis().SetTitleSize(0.07)
    draw_th1d(rat_sbe, rt.kFullCircle, rt.kBlue, 1.0, first=True,
              legend=leg3, text='(1/SBE:pp)/(1/SBE:Quenched)')
    rat_sbe.GetYaxis().SetRangeUser(0., 1.5)
    rat_sbe.Draw("PE")

    rat_5.Draw("PE same")
    rat_15.Draw("PE same")
    rat_50.Draw("PE same")
    leg_star3 = rt.TLegend(0.20, 0.70, .5, .95)
    leg_star3.AddEntry(rat_5, '(RooUnfold 5)/(Quench Truth)', 'p')
    leg_star3.AddEntry(rat_15, '(RooUnfold 15)/(Quench Truth)', 'p')
    leg_star3.AddEntry(rat_50, '(RooUnfold 50)/(Quench Truth)', 'p')
    # leg_star3.Draw()
    # leg3.Draw()

    c.SaveAs(f'{dir_name}/{outfile}')
    
def plot_RAA_set(outname, infiles, tags, colors, shapes, alphas):
    i = 0

    fout = rt.TFile.Open('__plot_RAA_set.root','recreate')
    hg_truth = None
    hgs = []
    ers = []
    for file, shape, color, alpha in zip(infiles, shapes, colors, alphas):
        i = i+1
        fin = rt.TFile.Open(file)
        _hg = fin.Get('rat_corr_AuAu_ruu4')
        _er = fin.Get('err_rat_bayes_AuAu')
        fout.cd()
        hg=_hg.Clone(f"hg_{i}")
        er=_er.Clone(f"er_{i}")
        # print(f'hg mean: {hg.GetMean()}')
        hg.SetMarkerStyle(shape)
        # print(f'color alpha {color} {alpha}')
        hg.SetMarkerColorAlpha(color, alpha+0.4)
        hg.SetLineColorAlpha(color, alpha+0.4)
        hg.SetMarkerStyle(shape)
        er.SetFillColorAlpha(color, alpha-0.4)
        hg.Write()
        er.Write()
        hgs.append(hg)
        ers.append(er)

        if not hg_truth:
            _hg_truth = fin.Get('rat_truth_AuAu')
            fout.cd()
            hg_truth = _hg_truth.Clone('hg_truth')
            hg_truth.SetMarkerStyle(rt.kOpenSquare)
            hg_truth.SetMarkerColorAlpha(rt.kRed, 1.)
            hg_truth.SetLineColorAlpha(rt.kRed, 1.)
            hg_truth.Write()
    fout.cd()

    canv = rt.TCanvas('c', 'c', 500, 400)
    canv.SetLeftMargin(0.20)
    canv.SetBottomMargin(0.20)
    canv.Draw()
    # leg = rt.TLegend(0.62, 0.7, 0.85, 0.87)

    hg_truth.SetTitle("")
    hg_truth.GetXaxis().SetRangeUser(16, 46)
    hg_truth.SetStats(0)
    hg_truth.GetYaxis().SetTitle("Ratio #frac{dN_{LeadJet}}{d#it{p}_{T,jet}^{truth}} to #it{pp}")
    hg_truth.GetYaxis().SetLabelSize(0.06)
    hg_truth.GetXaxis().SetTitle("#it{p}_{T,jet}^{truth} [GeV/#it{c}]")
    hg_truth.GetYaxis().SetRangeUser(0., 1.3)
    hg_truth.GetXaxis().SetRangeUser(16, 46)
    hg_truth.SetLineColor(kBlack)
    hg_truth.GetXaxis().SetTitleSize(0.07)
    hg_truth.GetXaxis().SetLabelSize(0.06)
    hg_truth.GetYaxis().SetTitleOffset(0.8)
    hg_truth.Draw("PE")
    leg.SetLineWidth(0)
    leg.AddEntry(hg_truth, r"Actual Ratio")

    for hg, er, tag in zip(hgs, ers, tags):
        hg.Draw("PE same")
        er.Draw("PE3 same")
        leg.AddEntry(hg, tag, 'p')

    leg.Draw()
    canv.SaveAs(outname) 
    canv.SaveAs(outname.replace('.png', '.pdf')) 


def plot_RAA_like(infile, pdf_name = None, rebins= None):
    if 'rhoAonly' in infile:
        tag = "NN_{AB}"
    elif 'reco_cpt' in infile:
        tag = 'NN_{pTcons}'
    elif 'reco_all' in infile:
        tag = 'NN_{RecoAll}'
    elif 'nconsts' in infile:
        tag = 'NN_{N_{cons}}'
    elif 'angularity' in infile:
        tag = 'NN_{Ang}'
    elif 'AB':
        tag = "AB Method"


    print(infile)
    fin = rt.TFile.Open(infile)
    if not fin:
        print(f'Cannot open file {infile}')
        sys.exit()
    else:
        print(f'Opened file {infile}')

    if pdf_name is None:
        pdf_name = infile.replace('.root', '.png')

    # 1. Create a TCanvas
    c = rt.TCanvas("c", "Canvas", 500, 500)
    c.SetLeftMargin(0.40)
    c.Draw()

    # 2. Divide the TCanvas into 2 TPads
    pad1 = rt.TPad("pad1", "Pad 1", 0, 0.6, 1, 1.0)
    pad1.SetBottomMargin(0)  # Upper and lower plot are joined
    pad1.SetLeftMargin(0.25)
    pad1.Draw()
    pad2 = rt.TPad("pad2", "Pad 2", 0, 0.06, 1, 0.6)
    pad2.SetTopMargin(0)
    pad2.SetBottomMargin(0.3)
    pad2.SetLeftMargin(0.25)
    pad2.Draw()

    # 3. In the top pad: Plot h2_pp and h2_quench
    pad1.cd()
    pad1.SetLogy()


    T_pp = fin.Get("truth_pp")
    T_AuAu = fin.Get("truth_AuAu")
    corr_AuAu_ruu4 = fin.Get("corr_AuAu_ruu_4")
    bayes_AuAu_error = fin.Get("bayes_AuAu_error")
    corr_AuAu_sbe = fin.Get("corr_AuAu_sbe")

    rat_corr_AuAu_sbe = fin.Get("rat_corr_AuAu_sbe")
    rat_corr_pp_sbe = fin.Get("rat_corr_pp_sbe")
    rat_corr_AuAu_sbe = fin.Get("rat_corr_AuAu_sbe")
    rat_corr_AuAu_ruu4 = fin.Get("rat_corr_AuAu_ruu4")
    err_rat_bayes_AuAu = fin.Get("err_rat_bayes_AuAu")
    rat_truth_AuAu = fin.Get("rat_truth_AuAu")

    if rebins:
        bins = np.array(rebins, dtype='d')
        rebin='rebin'
        T_pp = T_pp.Rebin(len(rebins)-1, 'T_pp_{rebin}', bins)
        T_AuAu = T_AuAu.Rebin(len(rebins)-1, f'T_AuAu_{rebin}', bins)
        corr_AuAu_sbe = corr_AuAu_sbe.Rebin(len(rebins)-1, f'corr_AuAu_sbe_{rebin}', bins)
        corr_AuAu_ruu4 = corr_AuAu_ruu4.Rebin(len(rebins)-1, f'corr_AuAu_ruu4_{rebin}', bins)
        bayes_AuAu_error = bayes_AuAu_error.Rebin(len(rebins)-1, f'bayes_AuAu_error_{rebin}', bins)

        T_pp.Scale(1./T_pp.GetBinWidth(1))
        T_AuAu.Scale(1./T_AuAu.GetBinWidth(1))
        corr_AuAu_sbe.Scale(1./corr_AuAu_sbe.GetBinWidth(1))
        corr_AuAu_ruu4.Scale(1./corr_AuAu_ruu4.GetBinWidth(1))
        bayes_AuAu_error.Scale(1./bayes_AuAu_error.GetBinWidth(1))


        
        rat_truth_AuAu = single_bin_mult(T_AuAu, T_pp, f'rat_truth_AuAu_{rebin}', div=1)
        rat_corr_pp_sbe = single_bin_mult(corr_AuAu_sbe, T_pp, f'rat_corr_pp_sbe_{rebin}', div=1)
        rat_corr_AuAu_sbe = single_bin_mult(corr_AuAu_sbe, T_pp, f'rat_corr_AuAu_sbe_{rebin}', div=1)
        rat_corr_AuAu_ruu4 = single_bin_mult(corr_AuAu_ruu4, T_pp, f'rat_corr_AuAu_ruu4_{rebin}', div=1)


        corr_AuAu_ruu2 = fin.Get("corr_AuAu_ruu_2")
        corr_AuAu_ruu6 = fin.Get("corr_AuAu_ruu_6")

        corr_AuAu_ruu2 = corr_AuAu_ruu2.Rebin(len(rebins)-1, f'corr_AuAu_ruu2_{rebin}', bins)
        corr_AuAu_ruu6 = corr_AuAu_ruu6.Rebin(len(rebins)-1, f'corr_AuAu_ruu6_{rebin}', bins)

        corr_AuAu_ruu2.Scale(1./corr_AuAu_ruu2.GetBinWidth(1))
        corr_AuAu_ruu6.Scale(1./corr_AuAu_ruu6.GetBinWidth(1))

        bayes_AuAu_error = make_err(corr_AuAu_ruu4, [corr_AuAu_ruu2, corr_AuAu_ruu6], rt.kGreen+2, 0.5)

        rat_ruu2 = single_bin_mult(corr_AuAu_ruu2, T_pp, f'rat_corr_AuAu_ruu2_{rebin}', div=1)
        rat_ruu6 = single_bin_mult(corr_AuAu_ruu6, T_pp, f'rat_corr_AuAu_ruu6_{rebin}', div=1)
        err_rat_bayes_AuAu = make_err(rat_corr_AuAu_ruu4, [rat_ruu2, rat_ruu6], rt.kGreen+2, 0.5)

        pdf_name = pdf_name.replace('.png', '_rebin.png')

    corr_AuAu_sbe.SetMarkerStyle(rt.kOpenTriangleUp)
    T_pp.SetMarkerStyle(rt.kFullCircle)
    T_pp.SetMarkerColor(rt.kBlack)
    T_pp.SetLineColor(rt.kBlack)
    T_pp.GetYaxis().SetTitle("#frac{dN_{LeadJet}}{d#it{p}_{T}}")
    T_pp.GetYaxis().SetLabelSize(0.10)
    T_pp.GetXaxis().SetRangeUser(corr_AuAu_sbe.GetXaxis().GetBinLowEdge(1), corr_AuAu_sbe.GetXaxis().GetBinUpEdge(corr_AuAu_sbe.GetNbinsX()))
    T_pp.GetYaxis().SetRangeUser(0.399e-8,6.0e-5)
    T_pp.GetXaxis().SetTitle("#it{p}_{T}")
    T_pp.GetXaxis().SetTitleOffset(1.4)
    T_pp.GetYaxis().SetTitleSize(0.10)
    T_pp.SetStats(0)
    T_pp.GetYaxis().SetTitleOffset(0.8)
    T_pp.SetMarkerColor(rt.kBlack)
    T_pp.SetLineColor(rt.kBlack)
    T_pp.SetTitle("")
    T_pp.Draw("PE1")
    # T_pp.Draw("E3 same")

    T_AuAu.SetMarkerStyle(rt.kOpenSquare)
    T_AuAu.SetMarkerColor(rt.kRed)
    T_AuAu.SetLineColor(rt.kRed)
    T_AuAu.Draw("PE SAME")

    corr_AuAu_sbe.SetMarkerColor(rt.kBlue)
    corr_AuAu_sbe.SetLineColor(rt.kBlue)
    corr_AuAu_sbe.SetMarkerSize(1.5)
    corr_AuAu_sbe.Draw("PE SAME")

    corr_AuAu_ruu4.SetMarkerStyle(rt.kOpenCircle)
    corr_AuAu_ruu4.SetMarkerSize(1.3)
    corr_AuAu_ruu4.SetMarkerColor(rt.kGreen+2)
    corr_AuAu_ruu4.SetLineColor(rt.kGreen+2)
    corr_AuAu_ruu4.Draw("PE SAME")

    bayes_AuAu_error.Draw("E3 same")

    leg = rt.TLegend(0.45, 0.55, 0.89, 0.88)
    leg.SetFillColorAlpha(rt.kWhite,0)
    leg.AddEntry(T_pp, "JETSCAPE #it{pp}", "p")
    leg.AddEntry(T_AuAu, "3.5 fm: Quenched", "p")
    leg.AddEntry(corr_AuAu_sbe, f"3.5 fm: {tag} Unfolded", "p")
    leg.AddEntry(corr_AuAu_ruu4, f"3.5 fm: {tag} 1-Bin Eff.", "p")
    leg.SetLineWidth(0)
    leg.Draw()

    pad2.cd()
    rat_truth_AuAu.SetMarkerStyle(rt.kOpenSquare)
    rat_truth_AuAu.SetMarkerColor(rt.kRed)
    rat_truth_AuAu.SetLineColor(rt.kRed)

    rat_truth_AuAu.GetYaxis().SetTitle("R_{AA}^{LeadJet}")
    rat_truth_AuAu.GetYaxis().SetLabelSize(0.08)
    rat_truth_AuAu.GetYaxis().SetTitleSize(0.10)
    rat_truth_AuAu.GetYaxis().SetTitleOffset(0.8)
    rat_truth_AuAu.GetYaxis().SetRangeUser(0., 0.99)

    rat_truth_AuAu.GetXaxis().SetTitle("#it{p}_{T} [GeV/#it{c}]")
    rat_truth_AuAu.GetXaxis().SetLabelSize(0.08)
    rat_truth_AuAu.GetXaxis().SetLabelOffset(0.01)
    rat_truth_AuAu.GetXaxis().SetRangeUser(corr_AuAu_sbe.GetXaxis().GetBinLowEdge(1), corr_AuAu_sbe.GetXaxis().GetBinUpEdge(corr_AuAu_sbe.GetNbinsX()))
    rat_truth_AuAu.GetXaxis().SetTitleSize(0.085)
    rat_truth_AuAu.GetYaxis().SetLabelOffset(0.02)
    rat_truth_AuAu.GetXaxis().SetTitleOffset(1.1)
    rat_truth_AuAu.SetTitle("")
    rat_truth_AuAu.SetStats(0)
    rat_truth_AuAu.Draw("PE")

    rat_corr_pp_sbe.SetMarkerStyle(rt.kOpenCircle)
    rat_corr_pp_sbe.SetMarkerSize(1.5)
    rat_corr_pp_sbe.SetMarkerColor(rt.kRed)
    # rat_corr_pp_sbe.Draw("PE same")

    rat_corr_AuAu_sbe.SetMarkerStyle(rt.kOpenTriangleUp)
    rat_corr_AuAu_sbe.SetMarkerColor(rt.kBlue)
    rat_corr_AuAu_sbe.SetLineColor(rt.kBlue)
    rat_corr_AuAu_sbe.Draw("PE SAME")

    rat_corr_AuAu_ruu4.SetMarkerStyle(rt.kOpenCircle)
    rat_corr_AuAu_ruu4.SetMarkerColor(rt.kGreen+2)
    rat_corr_AuAu_ruu4.SetLineColor(rt.kGreen+2)
    rat_corr_AuAu_ruu4.Draw("PE SAME")

    err_rat_bayes_AuAu.SetMarkerStyle(rt.kDot)
    err_rat_bayes_AuAu.SetMarkerColorAlpha(rt.kGreen+2, 0.5)
    err_rat_bayes_AuAu.SetFillColorAlpha(rt.kGreen+2, 0.5)
    err_rat_bayes_AuAu.Draw("E3 same")
    
    c.Update()
    
    c.SaveAs(pdf_name)
    c.SaveAs(pdf_name.replace('.png', '.pdf'))
    print('done')

    
    
