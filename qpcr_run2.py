import pandas as pd
import tkinter as tk
import numpy as np
from tkinter import filedialog
import opentrons.execute
from opentrons import simulate
import requests
import json
import operator

def uploadstringasfiletodropbox(name, string_data):

    name = "/" + name
    url = 'https://content.dropboxapi.com/2/files/upload'
    headers = {'Authorization': 'Bearer {{AUTHORIZATION_CODE}}',
               'Dropbox-API-Arg': '{\"path\": \"' + name + '\",\"mode\": \"add\",\"autorename\": true,\"mute\": false,\"strict_conflict\": false}',
               'Content-Type': 'application/octet-stream'
               }
    requests.post(url, headers=headers, data=str.encode(string_data))

def wellstring96fromindex(index) -> str:
    return chr(ord('A') + int(index/12)) + str(index%12 + 1)

root = tk.Tk()
root.withdraw()
qpcr_file_path = filedialog.askopenfilename()
qpcr_file_data = None

if qpcr_file_path.endswith(".xlsx"):
    qpcr_file_data = pd.read_excel(qpcr_file_path)
if qpcr_file_path.endswith(".csv"):
    qpcr_file_data = pd.read_csv(qpcr_file_path)

qpcr_file_data.dropna(inplace=True)

all_targets = np.sort(np.unique('/'.join(qpcr_file_data.iloc[:, 2]).split("/")))
well_pos_from_target = {}
col_num = 1
charcode = ord('A')

for target in all_targets:
    well_pos_from_target[target] = chr(charcode) + str(col_num)
    col_num = col_num + 1
    if col_num > 6:
        col_num = 1
        charcode = charcode + 1


pcr_reactions = pd.DataFrame(columns=['reaction_well', 'sample_well', 'sample_name', 'mastermix_well','mastermix_name'])
reaction_index = 1

for index, row in qpcr_file_data.iterrows():
    new_row = {}
    vals = row.values
    targets = str(vals[2]).split("/")
    for target in targets:
        new_row['sample_well'] = str(vals[0])
        new_row['sample_name'] = str(vals[1])
        new_row['mastermix_well'] = well_pos_from_target[target]
        new_row['mastermix_name'] = target

        pcr_reactions.loc[reaction_index] = new_row
        reaction_index = reaction_index + 1


pcr_reactions.sort_values(by=['mastermix_name'], inplace=True)


charcode = ord('H')
col_num = 1

for index, row in pcr_reactions.iterrows():

    row['reaction_well'] = chr(charcode) + str(col_num)

    if(chr(charcode) == "A"):
        charcode = ord('H')
        col_num = col_num + 1
    else:
        charcode = charcode - 1


pcr_output_sheet = pd.DataFrame(columns=['Row', 'Column', '*Target Name', '*Sample Name','*Biological Group'])

well_count = 0

while well_count < 96:

    ws = wellstring96fromindex(well_count)

    new_row = {'Row': ws[0], 'Column': ws[1:len(ws)]}

    reaction = pcr_reactions[pcr_reactions['reaction_well'] == ws].values
    if reaction.shape[0] > 0:
        new_row['*Target Name'] = reaction[0][4]
        new_row['*Sample Name'] = reaction[0][2]

    pcr_output_sheet.loc[well_count] = new_row
    well_count = well_count + 1

uploadstringasfiletodropbox("test_folder/pcr_output.csv", pcr_output_sheet.to_csv(index=False, header=True))


metadata={"apiLevel": "2.7"}

PROTOCOL = opentrons.execute.get_protocol_api('2.7')

P20_TIPS1 = PROTOCOL.load_labware('opentrons_96_tiprack_20ul', 4)
P20_TIPS2 = PROTOCOL.load_labware('opentrons_96_tiprack_20ul', 5)

P20_PIPETTE = PROTOCOL.load_instrument('p20_single_gen2', 'left', tip_racks=[P20_TIPS1, P20_TIPS2])
P20_PIPETTE.flow_rate.aspirate = 40
P20_PIPETTE.flow_rate.dispense = 40

SAMPLE_PLATE = PROTOCOL.load_labware('nest_96_wellplate_200ul_flat', 2)
QPCR_PLATE = PROTOCOL.load_labware('biorad_96_wellplate_200ul_pcr', 3)
MMX_TUBE_RACK = PROTOCOL.load_labware('opentrons_24_tuberack_nest_1.5ml_snapcap', 1)


#LOAD IN MASTER MIXES------------------------------------------
prev_mmx_well = None
for index, row in pcr_reactions.iterrows():

    if row['mastermix_well'] != prev_mmx_well:
        if P20_PIPETTE.has_tip:

            P20_PIPETTE.drop_tip()

        P20_PIPETTE.pick_up_tip()

    mmx_well = MMX_TUBE_RACK.wells_by_name()[row['mastermix_well']]
    rx_well = QPCR_PLATE.wells_by_name()[row['reaction_well']]

    P20_PIPETTE.transfer(10.5, mmx_well, rx_well, new_tip='Never')

    prev_mmx_well = row['mastermix_well']

P20_PIPETTE.drop_tip()


#--------------------------------------------------------------
#ADD SAMPLES---------------------------------------------------
#--------------------------------------------------------------
for index, row in pcr_reactions.iterrows():

    sample_well = SAMPLE_PLATE.wells_by_name()[row['sample_well']]
    reaction_well = QPCR_PLATE.wells_by_name()[row['reaction_well']]
    P20_PIPETTE.transfer(9.5, sample_well, reaction_well, mix_after=(2, 15))

for c in PROTOCOL.commands():
    print(c)