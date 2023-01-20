#!/usr/bin/env python3

# NOTE: This is mostly copied from:
#       <https://jerilkuriakose.medium.com/recover-corrupt-excel-file-xls-or-xlsx-using-python-2eea6bb07aae>
#       I've also made it a function.

#-------# Copied code starts here #--------------------------------------------#

# Changing the data types of all strings in the module at once
from __future__ import unicode_literals
# Used to save the file as excel workbook
# Need to install this library
from xlwt import Workbook
# Used to open to corrupt excel file
import io

def recover_corrupted_file(filename):
    """
    This function tries to recover a "corrupted dataframe". It's mostly copied
    from here:
    <https://jerilkuriakose.medium.com/recover-corrupt-excel-file-xls-or-xlsx-using-python-2eea6bb07aae>
    """

    # Opening the file using 'utf-8' encoding
    file1 = io.open(filename, "r", encoding="utf-8")
    data = file1.readlines()

    # Creating a workbook object
    xldoc = Workbook()
    # Adding a sheet to the workbook object
    sheet = xldoc.add_sheet("Sheet1", cell_overwrite_ok=True)
    # Iterating and saving the data to sheet
    for i, row in enumerate(data):
        # Two things are done here
        # Removeing the '\n' which comes while reading the file using io.open
        # Getting the values after splitting using '\t'
        for j, val in enumerate(row.replace('\n', '').split('\t')):
            sheet.write(i, j, val)

    # Generating a new filename:
    new_filename = filename.replace(".xls", "_RECOVERED.xls")
    # Saving the file as an excel file
    xldoc.save(new_filename)
    return new_filename

