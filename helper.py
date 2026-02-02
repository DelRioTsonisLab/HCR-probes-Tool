# Single cell: widgets + callback + display
import builtins
import re
import csv
import emoji
import os
import ipywidgets as widgets
from ipyfilechooser import FileChooser
from IPython.display import display, clear_output
from Bio import SeqIO
from makercb import maker  # Load IDT probe generator

# Append the local "bin" folder to PATH
current_dir = os.getcwd()
blast_bin = os.path.join(current_dir, "bin")
os.environ["PATH"] += os.pathsep + blast_bin

# Function
def load_idt():
    # get user namespace to set globals in the notebook
    user_ns = get_ipython().user_ns

    # Separator
    sep = widgets.HTML("<hr style='width:100%;margin:8px 0;'>")
    
    # Gene symbol input
    name_w      = widgets.Text(
        value="",
        placeholder='e.g OTX2',
        description="",
        style={'description_width': 'initial'}
    )
    amplifier_w = widgets.Dropdown(
        options=['B1', 'B2','B3','B4','B5','B7','B9','B10','B11','B13','B14','B15','B17'],
        description='Amplifier',
        style={'description_width': 'initial'}
    )
    db_w        = FileChooser(
        '.',
        #title='Select the FASTA file by clicking button',
        show_hidden=False
    )
    fullseq_w   = widgets.Textarea(
        value="",
        placeholder="""GAATGTACCCAAACCTGGGATAATTAAAGTCCTTGCAAGAGACTGGGGGAGAGATAGGTT
TCTTGAATTGTTTCTTTGTTTCTCACCCTTGGGCTTTGTTATCTGCATTATTTATTTAGC
CGAGGGGTTCTTTGTGTCTGCGAATGGCCCCAATCAAGTTTTGCTTGAGACAATTAGCTG
GTGCCCGGCCGGGAATCCTATGCAAATGCCCCGGGCGGCCGTACAATGCGGGCAGTTGAA
GGCACCTGGCGGCCCCCCCCTCCGCCCGACGGCACCCCCAGCACCGCCCGGGGCAGGCAG
GGCGCTCGGGCGGACCCCCCTCGCCCGGATTTCGGCCAACTCTCCTCCCCGCTCCAACTT
TAGCATGATGTCTTATCTTAAGCAACCACCTTACGCAGTCAATGGCCTGAGTCTCACCAC
CTCGGGCATGGATTTGTTGCATCCGTCCGTCGGTTATCCCGCCACCCCCCGAAAGCAGCG
GCGGGAGCGCACCACCTTCACCCGGGCGCAGCTGGACGTGCTGGAGGCCCTGTTCGCCAA
GACCCGCTACCCCGACATCTTCATGCGGGAGGAGGTGGCCTTGAAAATCAACCTGCCCGA
GTCCAGAGTGCAGGTGTGGTTCAAAAACCGCCGGGCCAAGTGCCGGCAGCAGCAGCAGCA
GCAGCAGAACGGGGGCCAGAACAAGGTGAGGCCGGCCAAAAAGAAGAACTCGCCGGCCCG
GGAAGTGAGCTCGGAGAGCGGGACCAGCGGGCAGTTCACACCCCCCTCCAGCACCTCGGT
CCCCACCATTTCCAGCAGCAGTGCCCCCGTGTCCATCTGGAGCCCGGCGTCCATCTCCCC
GCTCTCCGATCCCCTGTCCACTTCTTCCTCCTGCATGCAGAGGTCCTACCCCATGACCTA
CACCCAGGCATCAGGTTACAGCCAAGGATATGCCGGCTCGACCTCCTATTTCGGAGGGAT
GGACTGTGGATCTTATTTGACCCCTATGCACCACCAGCTCCCCGGACCGGGGGCCACCCT
GAGTCCCATGGGTGCCAATGCGGTCACCAGCCACCTCAACCAGTCTCCAGCCTCCCTCTC
CACCCAGGGCTATGGAGCCTCCAGTTTGGGCTTTAACTCGACCACCGATTGCTTGGATTA
TAAAGACCAAACCGCCTCCTGGAAGTTAAACTTCAATGCTGACTGCTTGGATTATAAAGA""",
        layout=widgets.Layout(width="45%", height="400px")
    )
    
    # Numeric inputs
    pause_w  = widgets.BoundedIntText(value=0, min=0, max=10000, step=1,
                                      description="Pairing Start", style={'description_width': 'initial'})
    polyAT_w = widgets.BoundedIntText(value=5, min=0, max=1000,
                                      description="polyAT", style={'description_width': 'initial'})
    polyCG_w = widgets.BoundedIntText(value=5, min=0, max=1000,
                                      description="polyCG", style={'description_width': 'initial'})
    numbr_w  = widgets.IntText(value=0, description="numbr")
    
    # Widget for intercepting maker's second prompt
    choice_w = widgets.BoundedIntText(
        value=1,
        min=1,
        description="Probe group #",
        style={'description_width': 'initial'}
    )
    
    # Boolean flags
    choose_w = widgets.Checkbox(
        value=True,
        description="Choose potential longest probe sets?",
        style={'description_width': 'initial'},
        layout=widgets.Layout(width='1000px')
    )
    BlastProbes_w = widgets.Checkbox(
        value=False,
        #disabled=True,
        description="Potential Probes",
        style={'description_width': 'initial'},
        layout=widgets.Layout(width='1000px')
    )
    dropout_w = widgets.Checkbox(
        value=False,
        description="Drop low quality Probes",
        style={'description_width': 'initial'},
        layout=widgets.Layout(width='1000px')
    )
    show_w = widgets.Checkbox(
        value=False,
        #disabled=True,
        description="Show detailed blast outputs",
        style={'description_width': 'initial'},
        layout=widgets.Layout(width='1000px')
    )
    report_w = widgets.Checkbox(
        value=True,
        #disabled=True,
        description="Display parameters chosen",
        style={'description_width': 'initial'},
        layout=widgets.Layout(width='1000px')
    )
    maxprobe_w = widgets.Checkbox(
        value=False,
        #disabled=True,
        description="Limit the number of max probes",
        style={'description_width': 'initial'},
        layout=widgets.Layout(width='1000px')
    )
    
    # Run button + output area
    run_button = widgets.Button(description="Store Inputs")
    output     = widgets.Output()
    
    
    
    
    ##########################
    # Option: manual or FASTA
    use_fasta = widgets.Checkbox(
        value=False, description='Select transcript from FASTA'
    )

    # Manual sequence entry
    seq_w = fullseq_w

    # FASTA file chooser and transcript controls
    fasta_w = FileChooser('.')
    get_tx = widgets.Button(description='Get transcripts')
    transcripts_info = widgets.HTML(value='')
    tx_dropdown = widgets.Dropdown(options=[], description='Transcript:')
    apply_tx = widgets.Button(description='Use transcript')
    preview = widgets.Output()

    # Hide controls initially
    tx_dropdown.layout.display = 'none'
    apply_tx.layout.display = 'none'
    transcripts_info.layout.display = 'none'
    
    # If use_fasta=False (manually-user)
    manual_box = widgets.VBox([
        widgets.HTML("<h4>Select the FASTA file to be used as the BLASTn subject:</h4>"),
        db_w
    ])
    # por defecto, si use_fasta=False (modo manual), lo mostramos:
    manual_box.layout.display = 'block'

    # Toggle manual vs FASTA mode
    def toggle_mode(change):
        if change['new']:
            seq_w.layout.display = 'none'
            fasta_box.layout.display = 'block'
            manual_box.layout.display   = 'none'
        else:
            fasta_box.layout.display = 'none'
            seq_w.layout.display = 'block'
            manual_box.layout.display   = 'block'

    use_fasta.observe(toggle_mode, names='value')

    # List transcripts matching gene in chosen FASTA
    def list_transcripts(b):
        transcripts_info.value = ''
        transcripts_info.layout.display = 'none'
        with preview:
            clear_output()
        # Hide dropdown and button
        tx_dropdown.layout.display = 'none'
        apply_tx.layout.display = 'none'
        tx_dropdown.options = []

        gene_input = name_w.value.strip()
        gene = gene_input.upper()
        filename = fasta_w.selected_filename
        directory = fasta_w.selected_path
        if not gene_input:
            with preview:
                print('❌ Provide gene name...')
            return
        if not filename:
            with preview:
                print('❌ Provide FASTA file...')
            return
        path = os.path.join(directory, filename)
        try:
            records = list(SeqIO.parse(path, 'fasta'))
        except Exception as e:
            with preview:
                print(f'❌ Error reading FASTA: {e}')
            return
        # Case-insensitive search in description
        matches = [rec.id for rec in records if gene in rec.description.upper()]
        if not matches:
            transcripts_info.value = '❌ No transcripts found.'
            transcripts_info.layout.display = 'block'
            return
        # Populate and show controls
        tx_dropdown.options = matches
        tx_dropdown.layout.display = 'block'
        apply_tx.layout.display = 'block'
        transcripts_info.value = f'Found {len(matches)} transcripts for {gene}.'
        transcripts_info.layout.display = 'block'

    # Load selected transcript into sequence box
    def load_transcript(b):
        with preview:
            clear_output()
        tid = tx_dropdown.value
        if not tid:
            return
        path = os.path.join(fasta_w.selected_path, fasta_w.selected_filename)
        for rec in SeqIO.parse(path, 'fasta'):
            if rec.id == tid:
                seq = str(rec.seq)
                # Remove spaces/newlines
                seq_w.value = re.sub(r'\s+', '', seq)
                with preview:
                    print(f'Loaded {rec.id} ({len(seq)} bp)')
                break

    get_tx.on_click(list_transcripts)
    apply_tx.on_click(load_transcript)

    # FASTA controls container (hidden by default)
    fasta_box = widgets.VBox([
        widgets.HTML('<b>Select cDNA FASTA file:</b>'),
        fasta_w,
        get_tx,
        transcripts_info,
        tx_dropdown,
        apply_tx,
        preview
    ])
    fasta_box.layout.display = 'none'


    def on_run_clicked(b):
        with output:
            clear_output()
            # 1) gene name
            if not name_w.value.strip():
                print("❌ Error: please enter a gene name.")
                return
            # 2)
            gene_input = name_w.value.strip()
            if not gene_input or not seq_w.value.strip():
                print('❌ Gene name or Sequence is missing...')
                return
            # Normalize gene name to uppercase
            user_ns['name']    = gene_input.upper()
            # Clean sequence
            user_ns['fullseq'] = re.sub(r'\s+', '', seq_w.value)
            user_ns['amplifier'] = amplifier_w.value

            # 3) FASTA file chosen?
            blast_flag = BlastProbes_w.value
            if blast_flag:
                # Use the same FASTA (fasta_w) that used in "Select transcripts"
                if use_fasta.value:
                    # fasta_w.selected es el FileChooser de tu FASTA de cDNA
                    user_ns['db'] = os.path.join(
                        fasta_w.selected_path,
                        fasta_w.selected_filename
                    )
                else:
                    # Use the FASTA (db_w) by selectting manually the transcript sequence
                    user_ns['db'] = db_w.selected
                # is there?
                # if not db_w.selected:
                #     print("❌ Error: FASTA invalid for BLASTn...")
                    # return
            else:
                # dont blast with BLASTn, FASTA is not enter
                user_ns['db'] = None

    
            # store every widget value into user_ns
            user_ns['name']       = name_w.value
            user_ns['fullseq']    = re.sub(r"\s+", "", fullseq_w.value)
            user_ns['amplifier']  = amplifier_w.value
            user_ns['pause']      = pause_w.value
            user_ns['choose']     = 'y' if choose_w.value else 'n'
            user_ns['polyAT']     = polyAT_w.value
            user_ns['polyCG']     = polyCG_w.value
            user_ns['BlastProbes']= 'y' if blast_flag else 'n' #'y' if BlastProbes_w.value else 'n'
            #user_ns['db']         = db_w.selected
            user_ns['dropout']    = 'y' if dropout_w.value else 'n'
            user_ns['show']       = 'y' if show_w.value else 'n'
            user_ns['report']     = 'y' if report_w.value else 'n'
            user_ns['maxprobe']   = 'y' if maxprobe_w.value else 'n'
            user_ns['numbr']      = numbr_w.value
            user_ns['choice']     = choice_w.value
    
            print(f"\u2714\ufe0f Inputs stored!")

    # bind the button
    run_button.on_click(on_run_clicked)

    # Display layout
    inputs = widgets.VBox([
        widgets.HTML("<h1>Custom Probe Design Parameters</h1>"),
        widgets.HTML("<h4>Gene Symbol:</h4>"),
        widgets.HBox([name_w]),
        sep,
        widgets.HTML("<h4>Enter your full sense-strand cDNA(5'→3'):</h4>"),
        seq_w,
        sep,
        use_fasta,
        fasta_box,
        sep,
        widgets.HTML("""<h4>Choose the hairpin you will use to amplify: <br>
        </h4><small>B1–B5 (Choi '14); B7–B17 (Wang '20)</small>"""),
        widgets.HBox([amplifier_w]),
        sep,
        widgets.HTML("<h4>How many bases from the 5' end of the Sense RNA before starting to hybridize?</h4>"),
        pause_w,
        widgets.HTML("<h4>Change the tolerated homopolymer lengths for poly (A & T) and poly (C & G):</h4>"),
        widgets.HBox([polyAT_w, polyCG_w]),
        sep,
        #widgets.HTML("<h4>Select the FASTA file to be used as the BLASTn subject:</h4>"),
        #db_w,
        manual_box,
        sep,
        widgets.VBox([
            widgets.HTML("""<h4>DEFAULTS PARAMETERS<br><br>
            Perform BLASTn on:</h4>"""),
            BlastProbes_w,
            widgets.HTML("<h4>Do you want to eliminate probes that appear in low quality BLAST outputs?</h4>"),
            dropout_w,
            widgets.HTML("<h4>Do you want to display detailed BLAST outputs?</h4>"),
            show_w,
            widgets.HTML("<h4>Do you want to display chosen parameters in output?</h4>"),
            report_w
        ]),
        sep,
        run_button,
        output
    ])
    
    display(inputs)
