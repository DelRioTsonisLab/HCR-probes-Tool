# Single cell: widgets + callback + display
import builtins
import re
import csv
import emoji  
import ipywidgets as widgets
from ipyfilechooser import FileChooser
from IPython.display import display, clear_output
from makercb import maker

def load_idt():
    # get user namespace to set globals in the notebook
    user_ns = get_ipython().user_ns
    
    # Separator
    sep = widgets.HTML("<hr style='width:100%;margin:8px 0;'>")
    
    # Text inputs
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
        title='Select your Fasta by clicking Button',
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
    pause_w  = widgets.BoundedIntText(value=0, min=0, max=200, step=1,
                                      description="pause", style={'description_width': 'initial'})
    polyAT_w = widgets.BoundedIntText(value=6, min=0, max=50,
                                      description="polyAT", style={'description_width': 'initial'})
    polyCG_w = widgets.BoundedIntText(value=6, min=0, max=50,
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
        value=True,
        #disabled=True,
        description="Blast Probes",
        style={'description_width': 'initial'},
        layout=widgets.Layout(width='1000px')
    )
    dropout_w = widgets.Checkbox(
        value=True,
        description="Drop low quality Probes",
        style={'description_width': 'initial'},
        layout=widgets.Layout(width='1000px')
    )
    show_w = widgets.Checkbox(
        value=True,
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
    
    def on_run_clicked(b):
        with output:
            clear_output()
            # 1) gene name
            if not name_w.value.strip():
                print("❌ Error: please enter a gene name.")
                return
            # 2) FASTA file chosen?
            if not db_w.selected:
                print("❌ Error: please select a FASTA file from the directory chooser.")
                return
    
            # store every widget value into user_ns
            user_ns['name']       = name_w.value
            user_ns['fullseq']    = re.sub(r"\s+", "", fullseq_w.value)
            user_ns['amplifier']  = amplifier_w.value
            user_ns['pause']      = pause_w.value
            user_ns['choose']     = 'y' if choose_w.value else 'n'
            user_ns['polyAT']     = polyAT_w.value
            user_ns['polyCG']     = polyCG_w.value
            user_ns['BlastProbes']= 'y' if BlastProbes_w.value else 'n'
            user_ns['db']         = str(db_w.selected)
            user_ns['dropout']    = 'y' if dropout_w.value else 'n'
            user_ns['show']       = 'y' if show_w.value else 'n'
            user_ns['report']     = 'y' if report_w.value else 'n'
            user_ns['maxprobe']   = 'y' if maxprobe_w.value else 'n'
            user_ns['numbr']      = numbr_w.value
            user_ns['choice']     = choice_w.value
    
            print(f"\u2714\ufe0f Inputs stored!")
    
    # bind the button
    run_button.on_click(on_run_clicked)
    
    # Layout & display
    inputs = widgets.VBox([
        widgets.HTML("<h1>HCR PROBE CUSTOM DESIGN</h1>"),
        widgets.HTML("<h4>Gene Name</h4>"),
        widgets.HBox([name_w]),
        sep,
        widgets.HTML("<h4>Enter your full sense-strand cDNA(5'→3') </h4>"),
        fullseq_w,
        sep,
        widgets.HTML("""<h4>Change the hairpin you will use to amplify with<br>
        (B1-B5 were used by Choi et al. 2014 and B7 to B17 were reported by Wang et al. BioRxiv 2020)</h4>"""),
        widgets.HBox([amplifier_w]),
        sep,
        widgets.HTML("<h4>How many bases from the 5' end of the Sense RNA before starting to hybridize?</h4>"),
        pause_w,
        widgets.HTML("<h4>Change the tolerated homopolymer lengths of (poly-A & poly-T) and/or (poly-C & poly-G)</h4>"),
        widgets.HBox([polyAT_w, polyCG_w]),
        sep,
        widgets.HTML("<h4>Specify the directory path of the FASTA file to be used as the BLASTn subject:</h4>"),
        db_w,
        sep,
        widgets.VBox([
            widgets.HTML("""<h4>Do you want to be able to select between potential longest probe sets?<br>
            (Unchecking defaults parameters)</h4>"""),
            choose_w,
            widgets.HTML("<h4>Do you want to eliminate probes that appear in low quality BLAST outputs?</h4>"),
            dropout_w,
            widgets.HTML("<h4>Toggle on/off whether a BLASTn search is performed on the potential probes or the original input cDNA:</h4>"),
            BlastProbes_w,
            #widgets.HTML("<h4>Do you want to limit the number of probes made</h4>"),
            #maxprobe_w,
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
