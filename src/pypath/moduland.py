#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright (c) 2014-2015 - EMBL-EBI
#
#  File author(s): Dénes Türei (denes@ebi.ac.uk)
#
#  Distributed under the GNU GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://www.ebi.ac.uk/~denes
#

#
# This module provides a python interface for the
# ModuLand community detection method family.
# Find details about the methods here:
# http://www.linkgroup.hu/modules.php
#

import os
import sys
import subprocess
import threading
import textwrap
import codecs
import datetime
import time
import random


class Moduland(object):
    def __init__(self,
                 network,
                 directed=None,
                 tmpdir='moduland_tmp',
                 weight=None,
                 name=None,
                 max_threads=None,
                 moduland_executables=None):
        '''
        This class is to run multiple ModuLand threads. 

        Parameters
        ----------
        pajekfile: str
            Path to the network in pajek format. You can get pajek 
            file from edgelist or igraph Graph object using the 
            update_network() method of the Moduland class.
        landscape: str
            The landscape method {linkland|nodeland|edgeweight|perturland}
        hill: str
            The gradient method {proportional|gradient|total}
        p: float
            The p parameter for PerturLand landscape method
        '''
        self.state = 'Initializing...'
        if os.name != 'posix':
            sys.stdout.write(
                'ModuLand executables can run only on POSIX'
                ' systems able to execute ELF32/System V/i386 format.\n')
            self.state = 'Error: non-posix system detected.'
            return None
        self.modexe = {
            'pajek_conv': 'pajek_conv',
            'linkland': 'linkland',
            'nodeland': 'nodeland',
            'edgeweight': 'edgeweight',
            'perturland': 'perturland',
            'total_lowmem': 'total_lowmem',
            'prop': 'prop',
            'cpxext_cat': 'cpxext_cat',
            'landscape_union': 'landscape_union',
            'mm.sh': 'mm.sh',
            'scor.sh': 'scor.sh',
            'dohist.sh': 'dohist.sh',
            'spearmerge.py': 'spearmerge.py',
            'modmerge': 'modmerge',
            'newlevel': 'newlevel',
            'projector': 'projector',
            'awk': 'awk',
            'ln': 'ln'
        }
        '''
        You are able to give paths or names of ModuLand progs,
        in case any of them is not in the known paths, or has 
        different name than default:
        '''
        if type(moduland_executables) is dict:
            for k, v in moduland_executables.iteritems():
                self.modexe[k] = v
        '''
        Checking if all the executables exist:
        '''
        for k, v in self.modexe.iteritems():
            if self.which(v) is None:
                self.state = 'Could not find executable: %s' % v
                sys.stdout.write(self.state + '\n')
                return None
        self.session = name if name is not None else self.gen_session_id()
        self.seq = 0
        self.threads = max_threads if max_threads is not None else self.cpu_count(
        )
        self.tmpdir = tmpdir
        if not os.path.exists(tmpdir):
            os.mkdir(tmpdir)
        self.fname_stems = []
        self.landscape_files = []
        self.hill_files = []
        self.processes = {}
        self.running = 0
        self.waiting = 0
        self.starting = 0
        self.free_threads = self.threads
        self.watching_thread = None
        '''
        If `network` is igraph object, and directed is None, 
        graph assumed to be directed, if the igraph network
        is directed. If `directed` is False, this overrides
        the igraph property, and the network will be imported
        as undirected.
        '''
        self.directed = directed
        self.update_network(network=network, weight=weight)
        self.landscape_commands = {
            'linkland': '%s' % (self.modexe['linkland']) + ' -0 %s %s',
            'nodeland': '%s' % (self.modexe['nodeland']) + ' -0 %s %s',
            'edgeweight':
            '%s' % (self.modexe['edgeweight']) + ' --in %s --out %s',
            'perturland': '%s' % (self.modexe['perturland']) + ' %s %f 0 %s 0'
        }
        self.gradient_commands = {
            'total': '%s' %
            (self.modexe['total_lowmem']) + ' --no-edge-belong -m 1000 %s',
            'gradient': '%s' %
            (self.modexe['prop']) + ' %s H --noEdgeBelong --gradientMethod',
            'proportional':
            '%s' % (self.modexe['prop']) + ' %s H --noEdgeBelong'
        }

    def clean_tmp(self, file_types=None, silent=False):
        for f in os.listdir(self.tmpdir):
            if f.startswith(self.session):
                if file_types is None or f[-3:] in file_types:
                    rm = os.path.join(self.tmpdir, f)
                    os.remove(rm)
                    if not silent:
                        sys.stdout.write('\t:: Removing %s\n' % rm)
                        sys.stdout.flush()

    def update_network(self, network, weight=None):
        self.network = network
        self.weight = weight
        if self.network.__class__.__name__ == 'list':
            self.pajek_from_edgelist()
        elif self.network.__class__.__name__ == 'Graph':
            self.pajek_from_igraph()
        else:
            self.console('Network should be an igraph Graph object '
                         'or a list of edges.')
        if self.pajek is not None:
            self.seq += 1
            self.fname_stems.append(
                os.path.join(self.tmpdir, self.session + '-%03d' % self.seq))
            self.pajekfile = self.fname_stems[-1] + '.net'
            self.write_pajek()
            self.get_cxg()

    def write_pajek(self):
        with codecs.open(self.pajekfile, encoding='utf-8', mode='w') as f:
            f.write(self.pajek)

    def pajek_from_edgelist(self):
        self.pajek = '*Vertices '
        self.vertices = []
        for e in network:
            if len(e) > 1:
                self.vertices.append(str(e[0]))
                self.vertices.append(str(e[1]))
            else:
                self.console(
                    'Edgelist contains element shorter than 2. '
                    'Each edge should have a source and a target node.')
                self.pajek = None
                return None
        self.vertices = list(set(self.vertices))
        self.vertices.sort()
        self.pajek += str(len(self.vertices))
        for i, v in enumerate(self.vertices):
            self.pajek += '%u %s\n' % (str(i + 1), v)
        if self.directed:
            self.pajek += '\n*Arcs\n'
        else:
            self.pajek += '\n*Edges\n'
        for e in network:
            self.edges.append((str(e[0]), str(e[1])))
            src = self.vertices.index(str(e[0])) + 1
            tgt = self.vertices.index(str(e[1])) + 1
            pajek_edge = [str(src), str(tgt)]
            if len(e) > 2:
                pajek_edge.append(str(e[2]))
                self.weight = True
            self.pajek += ' '.join(pajek_edge) + '\n'

    def pajek_from_igraph(self):
        self.pajek = '*Vertices %u\n' % self.network.vcount()
        self.vertices = []
        self.edges = []
        has_names = 'name' in self.network.vs.attributes()
        if self.directed is None:
            self.directed = self.network.is_directed()
        for v in self.network.vs:
            if has_names:
                self.pajek += str(v.index + 1)
            if has_names:
                self.pajek += ' %s' % v['name']
                self.vertices.append(v['name'])
            else:
                self.vertices.append(v.index)
            self.pajek += '\n'
        if self.directed:
            self.pajek += '*Arcs\n'
        else:
            self.pajek += '*Edges\n'
        for e in self.network.es:
            if has_names:
                self.edges.append((self.network.vs[e.source]['name'],
                                   self.network.vs[e.target]['name']))
            else:
                self.edges.append((e.source, e.target))
            pajek_edge = [str(e.source + 1), str(e.target + 1)]
            if self.weight is not None and self.weight in self.network.es.attributes(
            ):
                pajek_edge.append(str(e[self.weight]))
            self.pajek += ' '.join(pajek_edge) + '\n'

    def help(self):
        help_text = '''This module provides an interface to run the ModuLand community 
        detection methods family on igraph Graph objects. After instantiating the 
        class, you can start one or more ModuLand threads, monitor their run, and 
        get the results.
        More details, and citation: Kovács 2010, PLoS ONE 5(9):e12528.
        '''
        self.console(help_text)

    def get_cxg(self):
        self.cxg = self.pajekfile.replace('.net', '.cxg')
        cmd = ' '.join([self.modexe['pajek_conv'], self.pajekfile, self.cxg])
        self.console('Converting pajek file to cxg...')
        p = subprocess.Popen(cmd, shell=True)
        while True:
            time.sleep(0.5)
            returncode = p.poll()
            if returncode is not None:
                if os.path.exists(self.cxg):
                    self.console('Network has been converted to cxg format.')
                    self.console('Now ModuLand methods are available'
                                 'to run parallely.')
                    self.state = 'Ready to run landscape calculations.'
                else:
                    self.console('Failed to convert pajek file to cxg format.')
                    self.state = 'Failed to process the pajek file.'
                break

    def memberships(self, cxb):
        result = None
        cmd = ' '.join(
            [self.modexe['cpxext_cat'], cxb, '|', self.modexe['mm.sh']])
        proc = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        while True:
            time.sleep(0.5)
            state = proc.poll()
            if state is None:
                continue
            elif type(state) is int and state < 0:
                self.console('Failed to read memberships.\n\t%s' % cmd)
            else:
                result = proc.stdout.read()
                break
        return result

    def __str__(self):
        directed = 'directed' if self.directed else 'not directed'
        s = ' »» Moduland() instance\n'
        s += '    :: Network loaded: %u nodes, %u edges, %s\n' % (
            len(self.vertices), len(self.edges), directed)
        s += '    :: State: %s\n' % self.state
        s += '    :: Tasks: %u\n' % len(self.processes)
        for cmd, proc in self.processes.iteritems():
            this_cmd = proc['cmd']
            if len(this_cmd) > 80:
                this_cmd = '%s...\n               ...%s' % \
                    (this_cmd[:70], this_cmd[70:])
            s += '       [+] %s\n' % this_cmd
            if proc['state'] != 'waiting' and proc['state'] != 'starting':
                rtime = self.time_elapsed(
                    proc['started'], end=proc['finished'], string=True)
            else:
                rtime = '00:00:00'
            s += '           [state: %s :: running time: %s]\n' % \
                (proc['state'], rtime)
        return s

    def run_perturland(self, p=0.01, recalculate=False):
        if not self.weight:
            self.console('Unable to run perturland without edge weights.')
            return None
        lscape_cmd = {}
        cxl = self.fname_stems[-1] + '-perturland-%f.cxl' % p
        lscape_cmd['cmd'] = self.landscape_commands['perturland'] % (self.cxg,
                                                                     p, cxl)
        lscape_cmd['inputs'] = [self.cxg]
        lscape_cmd['outputs'] = [cxl]
        lscape_cmd['recalculate'] = recalculate
        self.run_anything(lscape_cmd)

    def run_landscape(self, landscape='linkland', p=0.01, recalculate=False):
        if landscape == 'perturland':
            self.run_perturland(p, recalculate=recalculate)
        else:
            lscape_cmd = {}
            cxl = self.fname_stems[-1] + '-' + landscape + '.cxl'
            lscape_cmd['cmd'] = self.landscape_commands[landscape] % (self.cxg,
                                                                      cxl)
            lscape_cmd['inputs'] = [self.cxg]
            lscape_cmd['outputs'] = [cxl]
            lscape_cmd['recalculate'] = recalculate
            self.run_anything(lscape_cmd)

    def run_hill(self, landscape_file, hill='proportional', recalculate=False):
        hill_cmd = {}
        cxb = landscape_file.replace('.cxl', '') + '-' + hill + '.cxb'
        hill_cmd['cmd'] = self.gradient_commands[hill] % \
            (' '.join([self.cxg, landscape_file, cxb]))
        hill_cmd['inputs'] = [landscape_file]
        hill_cmd['outputs'] = [cxb]
        hill_cmd['recalculate'] = recalculate
        self.run_anything(hill_cmd)

    def run_anything(self, param):
        if param['cmd'] in self.processes:
            if self.processes[param['cmd']]['process'].__class__.__name__ == 'Popen' \
                    and self.processes[param['cmd']]['state'] == 'running':
                self.console(
                    'Won\'t start command:\n\t%s\nThe same command already'
                    ' running, pid: %u' %
                    (param['cmd'],
                     self.processes[param['cmd']]['process'].pid))
                return None
            if self.processes[param['cmd']]['finished'] is not None \
                    and not param['recalculate']:
                sys.stdout.write(
                    'Command already run and finished:\n'
                    '\t%s\nTo run again, add `recalculate = True` parameter.\n'
                    % (param['cmd']))
                return None
        if self.check_files(param, silent=False):
            param['state'] = 'starting'
            self.starting += 1
        else:
            param['state'] = 'waiting'
            self.waiting += 1
        self.exec_cmd(param)

    def check_files(self, param, silent=True):
        msg = []
        missing = True
        overwrite = True
        for f in param['inputs']:
            if not os.path.exists(f):
                missing = False
                if not self.expect_output(f):
                    msg.append(' :: Input file %s missing.\n'
                               'Necessary input of scheduled command:\n%s\n'
                               'No task providing this file scheduled.\n' %
                               (f, param['cmd']))
        for f in param['outputs']:
            if os.path.exists(f) and not param['recalculate']:
                overwrite = False
                msg.append('Overwriting file %s' % f)
        if len(msg) > 0 and not silent:
            sys.stdout.write('\n'.join(msg))
            sys.stdout.flush()
        if param['recalculate']:
            return missing
        else:
            return missing and overwrite

    def expect_output(self, fname):
        '''
        Looking for scheduled tasks providing file `fname`.
        '''
        for proc in self.processes.values():
            if fname in proc['outputs']:
                return True
        return False

    def time_elapsed(self, started, end=None, string=False):
        end = datetime.datetime.now() if end is None else end
        diff = end - started
        hours, remainder = divmod(diff.seconds, 3600)
        minutes, seconds = divmod(remainder, 60)
        if string:
            if diff.days == 0:
                return '%02d:%02d:%02d' % \
                    (hours, minutes, seconds)
            else:
                return '%u days, %02d:%02d:%02d' % \
                    (diff.days, hours, minutes, seconds)
        else:
            return diff

    def exec_cmd(self, param):
        if self.watching_thread.__class__.__name__ != 'Thread' \
                or not self.watching_thread.isAlive():
            self.start_watching()
        '''
        Schedules an external process in `Moduland.processes`.
        '''
        self.processes[param['cmd']] = {
            'cmd': param['cmd'],
            'inputs': param['inputs'],
            'outputs': param['outputs'],
            'process': None,
            'state': param['state'],
            'started': None,
            'scheduled': datetime.datetime.now(),
            'finished': None,
            'recalculate': param['recalculate']
        }

    def start_process(self, cmd):
        '''
        Starts an external process.
        '''
        if cmd in self.processes and self.processes[cmd]['process'] is None:
            self.processes[cmd]['process'] = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                shell=True)
            self.processes[cmd]['state'] = 'running'
            self.processes[cmd]['started'] = datetime.datetime.now()

    def start_watching(self):
        self.wstop = False
        if self.watching_thread.__class__.__name__ != 'Thread' \
                or not self.watching_thread.isAlive():
            self.watching_thread = threading.Thread(target=self.watching)
            self.watching_thread.start()

    def stop_watching(self):
        self.wstop = True

    def watching(self):
        '''
        This is a scheduler to run the appropriate number of threads 
        with the proper order of execution of tasks.
        '''
        while True:
            time.sleep(0.5)
            for cmd, proc in self.processes.iteritems():
                if proc['finished'] is None:
                    if proc['state'] == 'waiting':
                        if self.check_files(proc):
                            self.processes[cmd]['state'] = 'starting'
                            self.waiting -= 1
                            self.starting += 1
                    if proc['state'] == 'starting':
                        if self.free_threads > 0:
                            self.start_process(cmd)
                            self.free_threads -= 1
                            self.running += 1
                            self.starting -= 1
                    if proc['state'] not in ['starting', 'waiting']:
                        state = proc['process'].poll()
                        if state is None:
                            self.processes[cmd]['state'] = 'running'
                            self.processes[cmd]['finished'] = None
                        else:
                            if self.processes[cmd]['finished'] is None:
                                if type(state) is int and state < 0:
                                    self.processes[cmd]['state'] = 'killed'
                                    self.running -= 1
                                elif sum([os.path.exists(x)
                                          for x in proc['outputs']]) < \
                                        len(proc['outputs']):
                                    self.processes[cmd]['state'] = 'failed'
                                    self.running -= 1
                                else:
                                    self.processes[cmd]['state'] = 'ready'
                                    self.running -= 1
                                self.processes[cmd]['returncode'] = state
                                self.processes[cmd]['stdout'] = proc[
                                    'process'].stdout
                                self.processes[cmd]['stderr'] = proc[
                                    'process'].stderr
                                self.processes[cmd]['finished'] = \
                                    datetime.datetime.now()
            self.free_threads = self.threads - self.running
            if self.running > 0 or self.waiting > 0:
                self.state = 'Working: \n\t:: running processes: %u\n\t '\
                    ':: ready to start: %u\n\t:: waiting: %u\n\t '\
                    ':: max threads: %u\n' % \
                    (self.running, self.starting, self.waiting, self.threads)
            else:
                self.state = 'No tasks, nothing to do.'
                self.wstop = True
            if self.wstop:
                break

    def cpu_count(self):
        """
        Detects the number of CPUs on a system. Cribbed from pp.
        Source:
        http://codeliberates.blogspot.co.uk/2008/05/detecting-cpuscores-in-python.html
        """
        # Linux, Unix and MacOS:
        if hasattr(os, "sysconf"):
            if os.sysconf_names.has_key("SC_NPROCESSORS_ONLN"):
                # Linux & Unix:
                ncpus = os.sysconf("SC_NPROCESSORS_ONLN")
                if isinstance(ncpus, int) and ncpus > 0:
                    return ncpus
            else:  # OSX:
                try:
                    return int(os.popen2("sysctl -n hw.ncpu")[1].read())
                except:
                    return 1
        # Windows:
        if os.environ.has_key("NUMBER_OF_PROCESSORS"):
            ncpus = int(os.environ["NUMBER_OF_PROCESSORS"])
            if ncpus > 0:
                return ncpus
        return 1  # Default

    def which(self, exe):
        '''
        Checks if executable is available.
        Source:
        http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
        '''

        def is_exe(fpath):
            return os.path.exists(fpath) and os.access(fpath, os.X_OK)

        def ext_candidates(fpath):
            yield fpath
            for ext in os.environ.get("PATHEXT", "").split(os.pathsep):
                yield fpath + ext

        fpath, fname = os.path.split(exe)
        if fpath:
            if is_exe(exe):
                return exe
        else:
            for path in os.environ["PATH"].split(os.pathsep):
                exe_file = os.path.join(path, exe)
                for candidate in ext_candidates(exe_file):
                    if is_exe(candidate):
                        return candidate
        return None

    # These duplicates are here to make this module independent
    # from other parts of `pypath`:

    def console(self, message):
        message = '\n\t'.join(textwrap.wrap(message, 50))
        sys.stdout.write(('\n\t' + message + '\n\n').ljust(50))
        sys.stdout.flush()

    def gen_session_id(self, length=5):
        abc = '0123456789abcdefghijklmnopqrstuvwxyz'
        return ''.join(random.choice(abc) for i in range(length))
