#!/home/davidstewart/penv-ML/bin/python3

import gzip
import ROOT
import sys
from os import path
from ROOT import std
from array import array
import time

def gunziptree(in_file):
    print(f'starting gunziptree {in_file}')
    start_time = time.time()

    if not path.isfile(in_file):
        print('File not found: {}'.format(in_file))
        sys.exit(1)


    # make the output file
    f_outname = in_file.replace('.dat.gz', '.root')

    f = ROOT.TFile(f_outname, 'recreate')
    tree = ROOT.TTree('T', 'tree')


    event = array('i',[-1])
    sigma = array('d',[-1.0])
    sigma_error = array('d',[-1.0])

    # jet hadrons
    # num = std.vector('int')()
    pid = std.vector('int')()
    pt = std.vector('float')()
    eta = std.vector('float')()
    phi = std.vector('float')()
    E = std.vector('float')()

    # energy loss show initiating parton ("ip")
    ip_pid = std.vector('int')()
    ip_pt = std.vector('float')()
    ip_eta = std.vector('float')()
    ip_phi = std.vector('float')()
    ip_E = std.vector('float')()

    # tree.Branch('num', num)
    tree.Branch('event', event,'event/I')
    tree.Branch('sigma', sigma,'sigma/D')
    tree.Branch('sigma_error', sigma_error,'sigma_error/D')
    tree.Branch('PID', pid)
    tree.Branch('pt', pt)
    tree.Branch('eta', eta)
    tree.Branch('phi', phi)
    tree.Branch('E', E)

    tree.Branch('ip_PID', ip_pid)
    tree.Branch('ip_pt', ip_pt)
    tree.Branch('ip_eta', ip_eta)
    tree.Branch('ip_phi', ip_phi)
    tree.Branch('ip_E', ip_E)

    # itest = 0
    read_next_as_initiating_shower = False
    is_in_final_state_hadrons = False
    i_hadron = -1

    for line in gzip.open(in_file,'rt'):
        # itest += 1
        # if itest % 10000 == 0:
            # print('Line: {}'.format(itest))
        # if itest > 200:
            # break

        if read_next_as_initiating_shower:
            read_next_as_initiating_shower = False
            arr = line.split()
            ip_pid.push_back(int(arr[1]))
            ip_pt.push_back(float(arr[3]))
            ip_eta.push_back(float(arr[4]))
            ip_phi.push_back(float(arr[5]))
            ip_E.push_back(float(arr[6]))
            continue

        if is_in_final_state_hadrons:
            arr = line.split()
            if len(arr) < 2 or arr[1] != 'H':
                is_in_final_state_hadrons = False
                i_hadron = -1

            else:
                if int(arr[0][1:-1]) != (i_hadron+1) or len(arr) != 13:
                    print('Error in reading final state hadron in line:')
                    print(line)
                    exit()
                i_hadron += 1

                _pt = float(arr[5])
                _eta = float(arr[6])
                if (_pt >= 0.2) and (abs(_eta) <=1.1):
                    pid.push_back(int(arr[3]))
                    pt.push_back(float(_pt))
                    eta.push_back(float(_eta))
                    phi.push_back(float(arr[7]))
                    E.push_back(float(arr[8]))
                    # num.push_back(pid.size()-1)
                    
                    # if (float(arr[8])>80):
                        # print(line)
                continue

        if line == '# Energy loss Shower Initating Parton: JetEnergyLoss\n':
            read_next_as_initiating_shower = True

        elif line == '# Final State Hadrons\n':
            is_in_final_state_hadrons = True

        elif line.endswith('Event\n'):
            _event = int(line.split()[0])
            if _event % 10000 == 0:
                print('Event: {}'.format(_event))
            if _event > 0:
                tree.Fill()

                pid.clear()
                pt.clear()
                eta.clear()
                phi.clear()
                E.clear()
                # num.clear()

                ip_pid.clear()
                ip_pt.clear()
                ip_eta.clear()
                ip_phi.clear()
                ip_E.clear()
            event[0] = _event

        elif line.startswith('# JetScape writersigmaGen'):
            sigma[0] = float(line.split()[-1])

        elif line.startswith('# JetScape writersigmaErr'):
            sigma_error[0] = float(line.split()[-1])

        elif line == '# Energy loss Shower Initating Parton: JetEnergyLoss':
            read_next_as_initiating_shower = True

    # if pid.size() > 0:
    tree.Fill()
    f.Write()

    print('Time taken: {:.2f} seconds'.format(time.time() - start_time))

if __name__ == '__main__':
    if len(sys.argv) > 1:
        in_file = sys.argv[1]
    else:
        in_file = 'test.dat.gz'

    gunziptree(in_file)
