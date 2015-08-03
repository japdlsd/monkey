import configparser
import sys
import signal
import re
import copy

CONFIG_PATH = 'config.ini'
# global flag for stopping program by SIGINT (CTRL+C on Linux)
isStopped = False 

# STATUSES {{{
OPEN = 'open'
OVERLAPPING = 'overlapping'
TOO_LONG = 'too_long'
STUCK = 'stuck' # not finished, but cannot continue
CLOSED = 'closed'
# }}}

def defaultConfig():
    config = configparser.ConfigParser()
    config['PATHS'] = {'database_path': 'database.fa',
                       'supercontigs_path': 'supercontigs.txt',
                       'supercontigs_output_path': 'supercontigs_output.txt'}
    config['SEARCHER'] = {}
    config['ENHANCER'] = {'max_contig_length': '80000',
                          'max_suffix_length': '50',
                          'min_suffix_length': '30',
                          'suffix_length_step': '5',
                          'successor_length': '8',
                          'definitive_successor_threshold': '0.95',
                          'definitive_successor_total_min': '40',
                          'branching_successor_threshold': '0.30',
                          'branching_successor_max_count': '3',
                          'branching_successor_total_min': '30',
                          'max_contig_amount': 5}
    config['SUPERCONFIGS'] = {'linebreak_at': '80',
                              'overlapping_min': '50',
                              'overlapping_max': '70'}
    return config

def verifyConfig(config):
    if int(config['SUPERCONFIGS']['linebreak_at']) <= 0: 
        print("!! linebreak_at should be positive number")

class Supercontigs:
    def __init__(self, config):
        self.config = config
        self.csp = re.compile('^\\s*$') # only whitespace characters
        self.op = re.compile('([ACGT]{'+ 
                        self.config['overlapping_min']
                        + ',' + 
                        self.config['overlapping_max']
                        + '})\#\\1') # overlapping abc#abc

    def read(self, filename):
        def is_comment(s):
            return len(s) > 0 and s[0] == '#'

        def is_status(s):
            return len(s) > 0 and s[0] == '@'
        
        def status(s):
            if is_status(s): return s[1:]
            raise Exception

        def is_contig_separator(s):
            nonlocal self
            return self.csp.match(s) is not None
       
        def is_supercontig_separator(s):
            return len(s) >= 5 and s[:5] == "-----"
        

        print("reading supercontigs from '{}'".format(filename))
        try:
            f = open(filename, 'r')
            
            self.array = []

            current_supercontig = []
            current_supercontig_status = OPEN
            current_contig = []
            for line in f.read().splitlines():
                # contigs separator
                if is_comment(line):
                    continue
                elif is_status(line):
                    current_supercontig_status = status(line)
                    if len(current_contig) > 0:
                        current_supercontig.append("".join(current_contig))
                        current_contig = []
                elif is_contig_separator(line):
                    if len(current_contig) > 0:
                        current_supercontig.append("".join(current_contig))
                        current_contig = []
                elif is_supercontig_separator(line):
                    if len(current_contig) > 0:
                        current_supercontig.append("".join(current_contig))
                        current_contig = []
                    if len(current_supercontig) > 0:
                        self.array.append({'status': current_supercontig_status, 
                                            'content': current_supercontig})
                        current_supercontig = []
                        current_supercontig_status = OPEN
                else:
                    current_contig.append(line)
            # finalizing
            if len(current_contig) > 0:
                current_supercontig.append("".join(current_contig))
            if len(current_supercontig) > 0:
                self.array.append({'status': current_supercontig_status, 
                                    'content': current_supercontig})

            f.close()
        except IOError as err:
            print(err)
            sys.exit(1)
        
        # add status for every contig and supercontig
        self.array = [{'status':sc['status'], 
                       'content':[ {'status':OPEN, 'content':contig} for contig in sc['content']]
                      } for sc in self.array]

    def write(self, filename):
        print("writing supercontigs to '{}'".format(filename))
        linebreaker = int(self.config['linebreak_at'])
        try:
            f = open(filename, 'w')
            for supercontig in self.array:
                print("-----", file=f)
                if supercontig['status'] != OPEN: print("@{}".format(supercontig['status']), file=f)
                for contig in supercontig['content']:
                    print(" ", file=f)
                    for i in range(0, len(contig['content']), linebreaker):
                        print(contig['content'][i:i+linebreaker], file=f)
            f.close()
        except IOError as err:
            print(err)
            sys.exit(1)

    def is_overlapping_contig(self, contig):
        for supercontig in self.array:
            if len(supercontig['content']) > 0:
                if self.op.match(contig['content'] + '#' + supercontig['content'][0]['content']) is not None:
                    return True
        return False

class Searcher:
    def __init__(self, config):
        self.config = config
    
    def read_database(self, filename):
        print("reading database from '{}' (could take a while)".format(filename))
        try:
            f = open(filename, 'r')
            self.database = []
            current_read = []
            for line in f.read().splitlines():
                if line[0] == '>':
                    self.database.append("".join(current_read))
                    current_read = []
                else:
                    current_read.append(line)
            self.database.append("".join(current_read))
            f.close()
        except IOError as err:
            print(err)
            sys.exit(1)

    def find_successors(self, suffix, successor_length):
        pattern = re.compile(suffix + '(.{' + str(successor_length) + '})', re.IGNORECASE)
        print(pattern)
        successors = {}
        for read in self.database:
            matched = pattern.findall(read)
            for succ in matched:
                if succ not in successors:
                    successors[succ] = 0
                successors[succ] += 1
        return successors


class Enhancer:
    def __init__(self, supercontigs, searcher, config):
        self.supercontigs = supercontigs
        self.searcher = searcher
        self.config = config
    
    def start(self):
        for i in range(len(self.supercontigs.array)):
            if isStopped: return
            supercontig = self.supercontigs.array[i]
            if supercontig['status'] != OPEN: continue
            print("supercontig number {}".format(i))
            
            stack = [x for x in supercontig['content'] if x['status'] == OPEN]
            supercontig['content'] = [x for x in supercontig['content'] if x['status'] != OPEN]
            
            total_contig_count = len(stack) + len(supercontig['content'])

            while len(stack) > 0:
                if isStopped: break
                contig = stack[-1]
                if self.supercontigs.is_overlapping_contig(contig):
                    contig['status'] = OVERLAPPING
                
                if len(contig['content']) > int(self.config['max_contig_length']):
                    contig['status'] = TOO_LONG

                if contig['status'] != OPEN:
                    stack.pop()
                    supercontig['content'].append(contig)
                    continue
                
                enhanced = False
                for suffix_length in range(int(self.config['max_suffix_length']), 
                                           int(self.config['min_suffix_length'])-1,
                                           -int(self.config['suffix_length_step'])):
                    if isStopped: break
                    successors = self.searcher.find_successors(contig['content'][-suffix_length:], 
                                                               int(self.config['successor_length']))
                    
                    successors = sorted([(y,x) for (x,y) in successors.items()])[::-1]
                    print(successors)
                    def definitive_successor(successors, threshold, total_minimum):
                        """successors are sorted!! (reversed)"""
                        if len(successors) == 0: return None
                        total = sum([a for (a,b) in successors])
                        if total < total_minimum: return None
                        if successors[0][0]/total >= threshold:
                            return successors[0][1]
                        else:
                            return None
                        
                    candidate = definitive_successor(successors,
                                                     float(self.config['definitive_successor_threshold']),
                                                     int(self.config['definitive_successor_total_min'])
                                                     )
                    if candidate != None:
                        print("definitive candidate found!")
                        contig['content'] += candidate
                        enhanced = True
                        break

                    def branching_successors(successors, threshold, total_minimum, max_count):
                        """successors are sorted!!! (reversed)"""
                        if len(successors) == 0: return None
                        total = sum([a for (a,b) in successors])
                        if total < total_minimum: return None
                        return [b for (a,b) in successors if a/total >= threshold][:max_count]
                
                    branching_candidates = branching_successors(
                                            successors,
                                            float(self.config['branching_successor_threshold']),
                                            int(self.config['branching_successor_total_min']),
                                            int(self.config['branching_successor_max_count'])
                                           )
                    if branching_candidates != None and len(branching_candidates) > 1:
                        print("branching candidates found!")
                        enhanced = True
                        config = stack.pop()
                        total_contig_count -= 1
                        for successor in branching_candidates:
                            if total_contig_count >= int(self.config['max_contig_amount']): break
                            new_contig = copy.deepcopy(contig)
                            new_contig['content'] += successor
                            stack.append(new_contig)
                            total_contig_count += 1

                    
                if not enhanced:
                    contig['status'] = STUCK
                    continue
            supercontig['content'] += stack

def main():
    print("Halo, mein lieber Freund :)")
    def signalHandler(signal, frame):
        global isStopped
        print("INTERRUPTED. Will end current operation and save the progress.")
        isStopped = True
    signal.signal(signal.SIGINT, signalHandler) 
    
    print("Importing settings from '{0}'".format(CONFIG_PATH))
    config = defaultConfig()
    config.read(CONFIG_PATH)
    
    verifyConfig(config)

    supercontigs = Supercontigs(config['SUPERCONFIGS'])
    supercontigs.read(config['PATHS']['supercontigs_path'])

    searcher = Searcher(config['SEARCHER'])
    searcher.read_database(config['PATHS']['database_path'])
    
    enhancer = Enhancer(supercontigs, searcher, config['ENHANCER'])

    enhancer.start()

    enhancer.supercontigs.write(config['PATHS']['supercontigs_output_path'])

if __name__ == "__main__":
    main()
