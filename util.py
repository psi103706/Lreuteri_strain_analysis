import subprocess

def run_process(cmd, stdout = subprocess.PIPE, stderr = subprocess.STDOUT, close_fds = True, shell = True, print_only = True):
    proc = subprocess.Popen(cmd, stdout = stdout, stderr = stderr, close_fds = close_fds, shell = shell)
    
    #only print output without parsing/piping lines
    if print_only:
        for line in iter(proc.stdout.readline, b''):
            print(line.decode('utf-8').rstrip())
        print()
        
    else:
        return proc