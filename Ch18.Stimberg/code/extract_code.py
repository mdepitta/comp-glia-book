from __future__ import print_function

import os, shutil
import glob

source_dir = '.'
target_dir = '../text/code'
example_dir = './clean_examples'

filenames = glob.glob(os.path.join(source_dir, 'example*.py'))
for filename in filenames:
    include = False
    ellipsis = False
    targets = []
    example_lines = []
    with open(filename, 'r') as f:
        for line in f:
            line = line[:-1]  # remove newline at end of file
            split_line = line.split('#', 1)
            if len(split_line) == 1:
                if include and not ellipsis:
                    targets[-1].append(line)
                example_lines.append(line)
            else:
                code, comment = split_line
                if 'INCLUDE BEGIN' in comment:
                    include = True
                    targets.append([])
                elif 'INCLUDE END' in comment:
                    include = False
                elif 'ELLIPSIS BEGIN' in comment:
                    ellipsis = True
                    # "code" should be only whitespace
                    targets[-1].append(code + '#  [...]')
                elif 'ELLIPSIS END' in comment:
                    ellipsis = False
                elif 'DELETE' in comment:
                    pass  # ignore the line completely
                else:
                    if include and not ellipsis:
                        targets[-1].append(line)
                    example_lines.append(line)
    print(filename, end=' ')

    if len(targets):  # anything?
        print('creating %d target files' % len(targets), end=' ')
        for i, target_lines in enumerate(targets):
            basename = os.path.basename(filename)
            without_ext, _ = os.path.splitext(basename)
            with open(os.path.join(target_dir, '%s_%d.py' % (without_ext, i+1)), 'w') as f:
                f.write('\n'.join(target_lines))
    with open(os.path.join(example_dir, os.path.basename(filename)), 'w') as f:
        f.write('\n'.join(example_lines))
    print()

# Copy over plot_utils and figures.mplstyle as well
additional_files = ['plot_utils.py', 'figures.mplstyle']
for filename in additional_files:
    shutil.copy(os.path.join(source_dir, filename), target_dir)
    shutil.copy(os.path.join(source_dir, filename), example_dir)
