'''
Created on 17.8.2015

@author: tohekorh
'''
import sys 

def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is one of "yes" or "no".
    """
    valid = {"yes":True,   "y":True,  "ye":True,
             "no":False,     "n":False}
    if default == None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "\
                             "(or 'y' or 'n').\n")


def write_line_own(wfile, line, key):
    log_f       =   open(wfile, key)
    log_f.write(line)            
    log_f.close()
    
def make_simul_param_file(fname, W, L, width, length_int, v, dy, T, \
                          dt, fric, thres_Z, interval, deltaY, theta, M, edge):

    log_f       =   open(fname, 'w')
    log_f.write('This is the params file. \n')
    
    log_f.write('width = %.6f [A] \n' %W)
    log_f.write('length = %.6f [A] \n' %L)
    log_f.write('width_i = %i \n' %width)
    log_f.write('length_i = %i \n' %length_int)
    log_f.write('Vel = %.6f [A/fs] \n' %v)
    log_f.write('Temp = %.1f [K] \n' %T)
    log_f.write('dt = %.4f [fs] \n' %dt)
    log_f.write('fric = %.4f \n' %fric)
    log_f.write('dy = %.6f [A] \n' %dy)
    log_f.write('thres_Z = %.4f [A] \n' %thres_Z)
    log_f.write('interval = %i \n' %interval)
    log_f.write('deltaY = %.6f [A] \n' %deltaY)
    log_f.write('theta = %.6f [A] \n' %theta)
    log_f.write('M = %i \n' %M)
    log_f.write('edge = %s \n' %edge)
    
    log_f.close()
    
def read_simul_params_file(fname):
    
    log_f =   open(fname, 'r')
    lines = log_f.readlines()
    
    width   =   find_between( lines[1],  ' = ', ' ' )
    length  =   find_between( lines[2],  ' = ', ' ' )
    width_i =   find_between( lines[3],  ' = ', ' ' )
    length_i=   find_between( lines[4],  ' = ', ' ' )
    vel     =   find_between( lines[5],  ' = ', ' ' )
    
    
    temp    =   find_between( lines[6],  ' = ', ' ' )
    dt      =   find_between( lines[7],  ' = ', ' ' )
    fric    =   find_between( lines[8],  ' = ', ' ' )
    dy      =   find_between( lines[9],  ' = ', ' ' )
    thresZ  =   find_between( lines[10], ' = ', ' ' )
    interval=   find_between( lines[11], ' = ', ' ' )
    deltaY  =   find_between( lines[12], ' = ', ' ' )
    theta   =   find_between( lines[13], ' = ', ' ' )
    M       =   find_between( lines[14], ' = ', ' ' )
    edge    =   find_between( lines[15], ' = ', ' ' )
    
    #width, length, width_i, length_i, v, T, dt, fric, dy, \
    #        thresZ, interval, deltaY theta, M, edge
    
    return float(width), float(length), int(width_i), int(length_i), float(vel), \
           float(temp), float(dt), float(fric), float(dy), \
            float(thresZ), int(interval), float(deltaY), float(theta), int(M), edge


def make_stick_simul_param_file(fname, width, length, length_dilde, T, \
                          dt, fric, interval, M, edge, stick):

    log_f       =   open(fname, 'w')
    log_f.write('This is the params file. \n')
    
    log_f.write('width_i = %i \n' %width)
    log_f.write('length = %.4f Angst \n' %length)
    log_f.write('lengthD = %.4f Angst  \n' %length_dilde)

    log_f.write('Temp = %.1f [K] \n' %T)
    log_f.write('dt = %.4f [fs] \n' %dt)
    log_f.write('fric = %.4f \n' %fric)
    log_f.write('interval = %i \n' %interval)
    log_f.write('M = %i \n' %M)
    log_f.write('edge = %s \n' %edge)
    log_f.write('stick = %s \n' %stick)
    
    
    log_f.close()
    
def read_stick_simul_params_file(fname):
    
    log_f =   open(fname, 'r')
    lines = log_f.readlines()
    
    width_i =   find_between( lines[1],  ' = ', ' ' )
    length  =   find_between( lines[2],  ' = ', ' ' )
    length_d=   find_between( lines[3],  ' = ', ' ' )
    temp    =   find_between( lines[4],  ' = ', ' ' )
    dt      =   find_between( lines[5],  ' = ', ' ' )
    fric    =   find_between( lines[6],  ' = ', ' ' )
    interval=   find_between( lines[7], ' = ', ' ' )
    M       =   find_between( lines[8], ' = ', ' ' )
    edge    =   find_between( lines[9], ' = ', ' ' )
    stick   =   find_between( lines[10], ' = ', ' ' )
    
    
    #width, length, width_i, length_i, v, T, dt, fric, dy, \
    #        thresZ, interval, deltaY theta, M, edge
    
    return int(width_i), float(length), float(length_d), float(temp), \
        float(dt), float(fric), int(interval), int(M), edge, stick in ['True', 'true', '1']
    
'''
def make_stick_simul_param_file(fname, width, length, length_dilde, T, \
                          dt, fric, interval, M, edge, stick):

    log_f       =   open(fname, 'w')
    log_f.write('This is the params file. \n')
    
    log_f.write('width_i = %i \n' %width)
    log_f.write('length = %.4f Angst \n' %length)
    log_f.write('lengthD = %.4f Angst  \n' %length_dilde)

    log_f.write('Temp = %.1f [K] \n' %T)
    log_f.write('dt = %.4f [fs] \n' %dt)
    log_f.write('fric = %.4f \n' %fric)
    log_f.write('interval = %i \n' %interval)
    log_f.write('M = %i \n' %M)
    log_f.write('edge = %s \n' %edge)
    log_f.write('stick = %s \n' %stick)
    
    
    log_f.close()
    
def read_stick_simul_params_file(fname):
    
    log_f =   open(fname, 'r')
    lines = log_f.readlines()
    
    width_i =   find_between( lines[3],  ' = ', ' ' )
    length  =   find_between( lines[4],  ' = ', ' ' )
    length_d=   find_between( lines[4],  ' = ', ' ' )
    temp    =   find_between( lines[6],  ' = ', ' ' )
    dt      =   find_between( lines[7],  ' = ', ' ' )
    fric    =   find_between( lines[8],  ' = ', ' ' )
    interval=   find_between( lines[11], ' = ', ' ' )
    M       =   find_between( lines[14], ' = ', ' ' )
    edge    =   find_between( lines[15], ' = ', ' ' )
    stick   =   find_between( lines[15], ' = ', ' ' )
    
    #width, length, width_i, length_i, v, T, dt, fric, dy, \
    #        thresZ, interval, deltaY theta, M, edge
    
    return int(width_i), float(length), float(length_d), float(temp), \
        float(dt), float(fric), int(interval), int(M), edge, stick in ['True', 'true', '1']
'''    
    
def find_between( s, first, last ):
    try:
        start = s.index( first ) + len( first )
        end = s.index( last, start )
        return s[start:end]
    except ValueError:
        return ""

def append_files(file1, file2, fname):
    
    fnew=   open(fname, 'w')
    f1  =   open(file1, 'r')
    f2  =   open(file2, 'r')
    
    lines = f1.readlines()
    for i in range(0, len(lines)):
        line = lines[i]
        fnew.write(line)
    
    lines = f2.readlines()
    for i in range(0, len(lines)):
        line = lines[i]
        fnew.write(line)
    
    fnew.close()
    f1.close()
    f2.close()
    