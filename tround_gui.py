"""This provides a GUI interface for tround.py."""

import platform
import ctypes
import os
import pickle
import io
from contextlib import redirect_stdout

import pyperclip
import PySimpleGUI as sg

import tround


def set_gui_colors():
    """Set the color scheme."""
    # this scheme is similar to Windows default colors
    sg.SetOptions(background_color='#F0F0F0',
                  text_element_background_color='#F0F0F0',
                  element_background_color='#F0F0F0',
                  text_color='#000000',
                  input_elements_background_color='#FFFFFF',
                  button_color=('#000000', '#E1E1E1'),
                  font=('Segoe UI', 9))
    print('Windows release is %s.' % (platform.release()))
    if int(platform.release()) >= 8:
        print('Registering DPI awareness.')
        ctypes.windll.shcore.SetProcessDpiAwareness(True)


def get_settings_filename():
    """Return the settings path and filename."""
    filename = os.path.join(os.environ['APPDATA'],
                            'tround_gui',
                            'config.ini')
    return filename


def load_values(window):
    """Load window fields from the settings file."""
    filename = get_settings_filename()
    # if file doesn't exist, don't load it
    if not os.path.isfile(filename):
        return
    # load it into a dictionary
    values = pickle.load(open(filename, 'rb'))
    print(values)
    for key, value in values.items():
        window.FindElement(key).Update(value)
    print('Loading window fields from disk')


def save_values(values):
    """Save window fields to the settings file."""
    filename = get_settings_filename()
    if not os.path.exists(os.path.dirname(filename)):
        os.makedirs(os.path.dirname(filename))
    print(values)
    print('Saving window fields to disk')
    pickle.dump(values, open(filename, 'wb'))


def process_board(window):
    """Perform the operations on the given file."""
    # set options based on checkboxes
    option = window.FindElement('round_corners').Get() == 1
    tround.rounded_corners = option
    option = window.FindElement('round_junctions').Get() == 1
    tround.rounded_junctions = option
    option = window.FindElement('teardrop_vias').Get() == 1
    tround.create_teardrops_on_vias = option
    option = window.FindElement('teardrop_pths').Get() == 1
    tround.create_teardrops_on_pths = option
    # read board filename
    filename = window.FindElement('filename').Get()
    layout = [[sg.Text('Performing operations on board file...')],
              [sg.Multiline(key='log', size=(80, 20))],
              [sg.Button('Copy command to clipboard', key='clipboard'),
               sg.OK(key='ok')]]
    status = sg.Window('Tround GUI').Layout(layout)
    status.Read(timeout=0)
    status.FindElement('log').Update('Reading board file')
    status.FindElement('ok').SetFocus()
    status.Read(timeout=0)
    status.TKroot.grab_set()
    # temp = builtins.print
    # builtins.print = sg.EasyPrint
    f = io.StringIO()
    with redirect_stdout(f):
        board = tround.Board(filename)
        status.FindElement('log').Update(f.getvalue())
        tround.round_signals(board)
        status.FindElement('log').Update(f.getvalue())
        tround.teardrop_board_vias(board)
        status.FindElement('log').Update(f.getvalue())
        board.backup_file()
        status.FindElement('log').Update(f.getvalue())
        script_filename = board.generate_script()
        status.FindElement('log').Update(f.getvalue())
    status.FindElement('log').TKText.see('end')
    # sg.EasyPrint(f.getvalue())
    # builtins.print = temp
    event, values = status.Read()
    status.TKroot.grab_release()
    status.Close()
    if event == 'clipboard':
        pyperclip.copy('script %s;' % script_filename)
        sg.Popup('Command was copied to the clipboard.',
                 'To run, open the board in Eagle and run the command.')


def run_gui():
    """Run the Tround GUI."""
    set_gui_colors()
    options_frame = [
        [sg.Checkbox('Round corners of traces',
                     default=True,
                     key='round_corners')],
        [sg.Checkbox('Round corners of junctions',
                     default=True,
                     key='round_junctions')],
        [sg.Checkbox('Create teardrops on vias',
                     default=True,
                     key='teardrop_vias')],
        [sg.Checkbox('Create teardrops on component PTHs',
                     default=True,
                     key='teardrop_pths')]]
    layout = [[sg.Text('Eagle board file'),
               sg.InputText(key='filename',
                            do_not_clear=True,
                            enable_events=True),
               sg.FileBrowse(file_types=(("Eagle board files", "*.brd"),
                                         ("All types", "*.*")))],
              [sg.Frame('Options', options_frame)],
              [sg.Button('Generate Script', key='ok'),
               sg.Button('Exit', key='cancel')]]
    window = sg.Window('Tround GUI').Layout(layout)
    window.Read(timeout=0)
    # load previous options
    load_values(window)
    # set focus to OK button
    window.FindElement('ok').SetFocus()
    # move cursor to end of filename entry and make it visible
    window.FindElement('filename').TKEntry.icursor('end')
    window.FindElement('filename').TKEntry.xview('end')
    while True:
        # set focus to OK button
        event, values = window.Read()
        # quit form on request
        if event is None or event == 'cancel':
            break
        # save values
        save_values(values)
        if event == 'filename':
            if window.FindElementWithFocus() is None:
                print('Reclaiming focus')
                window.FindElement('ok').SetFocus()
                window.FindElement('filename').TKEntry.selection_clear()
                # move cursor to end of filename entry and make it visible
                window.FindElement('filename').TKEntry.icursor('end')
                window.FindElement('filename').TKEntry.xview('end')
            continue
        elif event == 'ok':
            process_board(window)
            continue
        else:
            print('WARNING: unknown event: %s' % event)


# execute as a script
if __name__ == "__main__":
    run_gui()
