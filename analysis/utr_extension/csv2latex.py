#!/usr/bin/env python
#coding=utf-8
########################################
# only get the debug function if run from Ipython #
def run_from_ipython():
    try:
        __IPYTHON__
        return True
    except NameError:
        return False

if run_from_ipython():
    from IPython.Debugger import Tracer
    debug = Tracer()
else:
    def debug():
        print ('Warning: debugging mark present.')
        pass
########################################

import django
from django.template import Template, Context
import csv

if __name__ == "__main__":

    # This line is required for Django configuration
    try:
        django.conf.settings.configure()
    except RuntimeError:
        pass

    # Open and read CSV file
    fid = open("3utr_extension.tsf")
    reader = csv.reader(fid, delimiter='\t')

    # Open and read template
    with open("template.tex") as f:
        t = Template(f.read())

    # Define context with the table data
    head = reader.next()
    c = Context({"head": head, "table": reader})

    # Render template
    output = t.render(c)

    fid.close()

    # Write the output to a file
    with open("table.tex", 'w') as out_f:
        out_f.write(output)
