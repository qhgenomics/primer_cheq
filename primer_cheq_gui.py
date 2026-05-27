import sys, os
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QLabel, QLineEdit, QPushButton, QFileDialog, QCheckBox, QGridLayout, QMessageBox
from PyQt5.QtCore import QThread, pyqtSignal, QTimer
import subprocess
from primer_cheq import get_primer_sequences, get_primer_table, get_db_glob, get_db_fasta, get_db_fastas, get_db_folder, download_bac, download_virus, align_primers

class BlastWorker(QThread):
    finished = pyqtSignal()

    def __init__(self, primer_file, primer_dict, primer_info, working_directory, prefix, sassy_loc, max_indel, max_mismatch, indel_mult, threads, max_primer_dist):
        super().__init__()
        self.primer_file = primer_file
        self.primer_dict = primer_dict
        self.primer_info = primer_info
        self.working_directory = working_directory
        self.prefix = prefix
        self.sassy_loc = sassy_loc
        self.max_indel = max_indel
        self.max_mismatch = max_mismatch
        self.indel_mult = indel_mult
        self.threads = threads
        self.max_primer_dist = max_primer_dist

    def run(self):
        align_primers(self.primer_dict, self.primer_file, self.working_directory, self.prefix, self.max_indel, self.max_mismatch, self.indel_mult, self.threads, self.sassy_loc)
        self.finished.emit()

class PrimerCheqGUI(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        layout = QGridLayout()
        self.output_label = QLabel('Input options (required - use either Primer file OR Primer table)')
        self.output_label.setStyleSheet("background-color: lightyellow; font-weight: bold")
        layout.addWidget(self.output_label, 0, 0, 1, 4)
        # Primer file
        self.primer_label = QLabel('Primer FASTA file:')
        self.primer_input = QLineEdit(self)
        self.primer_button = QPushButton('Browse', self)
        self.primer_button.clicked.connect(self.browse_primer)
        layout.addWidget(self.primer_label, 1, 0)
        layout.addWidget(self.primer_input, 1, 1, 1, 2)
        layout.addWidget(self.primer_button, 1, 3)
        
        # Primer table
        self.primer_table_label = QLabel('Primer table file (TSV):')
        self.primer_table_input = QLineEdit(self)
        self.primer_table_button = QPushButton('Browse', self)
        self.primer_table_button.clicked.connect(self.browse_primer_table)
        layout.addWidget(self.primer_table_label, 2, 0)
        layout.addWidget(self.primer_table_input, 2, 1, 1, 2)
        layout.addWidget(self.primer_table_button, 2, 3)

        self.sassy_label = QLabel('Sassy location:')
        self.sassy_input = QLineEdit('sassy', self)
        self.sassy_button = QPushButton('Browse', self)
        self.sassy_button.clicked.connect(self.browse_sassy)
        layout.addWidget(self.sassy_label, 3, 0)
        layout.addWidget(self.sassy_input, 3, 1, 1, 2)
        layout.addWidget(self.sassy_button, 3, 3)

        self.output_label = QLabel('Output options (required)')
        self.output_label.setStyleSheet("background-color: lightyellow; font-weight: bold")
        layout.addWidget(self.output_label, 4, 0, 1, 4)
        # Prefix
        self.prefix_label = QLabel('Prefix:')
        self.prefix_input = QLineEdit(self)
        layout.addWidget(self.prefix_label, 5, 0)
        layout.addWidget(self.prefix_input, 5, 1, 1, 3)

        # Working directory
        self.working_dir_label = QLabel('Working directory:')
        self.working_dir_input = QLineEdit(self)
        self.working_dir_button = QPushButton('Browse', self)
        self.working_dir_button.clicked.connect(self.browse_working_dir)
        layout.addWidget(self.working_dir_label, 6, 0)
        layout.addWidget(self.working_dir_input, 6, 1, 1, 2)
        layout.addWidget(self.working_dir_button, 6, 3)

        self.output_label = QLabel('Genome database options (pick 1 or more)')
        self.output_label.setStyleSheet("background-color: lightyellow; font-weight: bold")
        layout.addWidget(self.output_label, 7, 0, 1, 4)
        # NCBI Virus
        self.ncbi_virus_label = QLabel('NCBI Virus Taxid (will download fastas):')
        self.ncbi_virus_input = QLineEdit(self)
        layout.addWidget(self.ncbi_virus_label, 8, 0)
        layout.addWidget(self.ncbi_virus_input, 8, 1, 1, 3)

        # NCBI Bacteria
        self.ncbi_bacteria_label = QLabel('NCBI Bacteria Taxid (will download fastas):')
        self.ncbi_bacteria_input = QLineEdit(self)
        layout.addWidget(self.ncbi_bacteria_label, 9, 0)
        layout.addWidget(self.ncbi_bacteria_input, 9, 1, 1, 3)

        # Year filter
        self.year_label = QLabel('Year filter (for NCBI downloads):')
        self.year_input = QLineEdit('all', self)
        layout.addWidget(self.year_label, 10, 0)
        layout.addWidget(self.year_input, 10, 1, 1, 3)

        # Directory DB
        self.directory_db_label = QLabel('Directory with multiple fasta files (1 file per genome):')
        self.directory_db_input = QLineEdit(self)
        self.directory_db_button = QPushButton('Browse', self)
        self.directory_db_button.clicked.connect(self.browse_directory_db)
        layout.addWidget(self.directory_db_label, 11, 0)
        layout.addWidget(self.directory_db_input, 11, 1, 1, 2)
        layout.addWidget(self.directory_db_button, 11, 3)

        # Fasta DB
        self.fasta_db_label = QLabel('Fasta file with multiple genomes:')
        self.fasta_db_input = QLineEdit(self)
        self.fasta_db_button = QPushButton('Browse', self)
        self.fasta_db_button.clicked.connect(self.browse_fasta_db)
        layout.addWidget(self.fasta_db_label, 12, 0)
        layout.addWidget(self.fasta_db_input, 12, 1, 1, 2)
        layout.addWidget(self.fasta_db_button, 12, 3)

        # Glob DB
        self.glob_db_label = QLabel('Path with wildcard (*) that resolves to one or more genomes:')
        self.glob_db_input = QLineEdit(self)
        layout.addWidget(self.glob_db_label, 13, 0)
        layout.addWidget(self.glob_db_input, 13, 1, 1, 3)

        self.advanced_label = QLabel('Advanced options (optional)')
        self.advanced_label.setStyleSheet("background-color: lightcyan; font-weight: bold")
        layout.addWidget(self.advanced_label, 14, 0, 1, 4)

        # Max indel
        self.max_indel_label = QLabel('Max indels:')
        self.max_indel_input = QLineEdit('2', self)
        layout.addWidget(self.max_indel_label, 15, 0)
        layout.addWidget(self.max_indel_input, 15, 1)

        # Max mismatch
        self.max_mismatch_label = QLabel('Max mismatches:')
        self.max_mismatch_input = QLineEdit('4', self)
        layout.addWidget(self.max_mismatch_label, 15, 2)
        layout.addWidget(self.max_mismatch_input, 15, 3)

        # Indel multiplier
        self.indel_mult_label = QLabel('Indel multiplier:')
        self.indel_mult_input = QLineEdit('2', self)
        layout.addWidget(self.indel_mult_label, 16, 0)
        layout.addWidget(self.indel_mult_input, 16, 1)

        # Threads
        self.threads_label = QLabel('Threads:')
        self.threads_input = QLineEdit('1', self)
        layout.addWidget(self.threads_label, 16, 2)
        layout.addWidget(self.threads_input, 16, 3)

        # Max primer distance
        self.max_primer_dist_label = QLabel('Max primer distance:')
        self.max_primer_dist_input = QLineEdit('5000', self)
        layout.addWidget(self.max_primer_dist_label, 17, 0)
        layout.addWidget(self.max_primer_dist_input, 17, 1)

        # Dataset location
        self.dataset_loc_label = QLabel('Dataset location:')
        self.dataset_loc_input = QLineEdit('datasets', self)
        self.dataset_loc_button = QPushButton('Browse', self)
        self.dataset_loc_button.clicked.connect(self.browse_dataset_loc)
        layout.addWidget(self.dataset_loc_label, 17, 2)
        layout.addWidget(self.dataset_loc_input, 17, 3)

        # Only product forming
        self.only_product_forming_checkbox = QCheckBox('Only report product-forming primers', self)
        layout.addWidget(self.only_product_forming_checkbox, 18, 0, 1, 4)

        # Run button
        self.run_button = QPushButton('Run', self)
        self.run_button.clicked.connect(self.run_script)
        layout.addWidget(self.run_button, 19, 0, 1, 4)

        self.setLayout(layout)
        self.setWindowTitle('Primer Cheq GUI')
        self.show()

    def browse_primer(self):
        fname = QFileDialog.getOpenFileName(self, 'Open file', '', 'FASTA files (*.fasta *.fa)')
        if fname[0]:
            self.primer_input.setText(fname[0])

    def browse_primer_table(self):
        fname = QFileDialog.getOpenFileName(self, 'Open file', '', 'TSV files (*.tsv *.txt)')
        if fname[0]:
            self.primer_table_input.setText(fname[0])

    def browse_sassy(self):
        fname = QFileDialog.getOpenFileName(self, 'Open file', '', 'Sassy binary (*.exe, *)')
        if fname[0]:
            self.sassy_input.setText(fname[0])

    def browse_working_dir(self):
        dname = QFileDialog.getExistingDirectory(self, 'Select directory')
        if dname:
            self.working_dir_input.setText(dname)

    def browse_directory_db(self):
        dname = QFileDialog.getExistingDirectory(self, 'Select directory')
        if dname:
            self.directory_db_input.setText(dname)

    def browse_fasta_db(self):
        fname = QFileDialog.getOpenFileName(self, 'Open file', '', 'FASTA files (*.fasta *.fa)')
        if fname[0]:
            self.fasta_db_input.setText(fname[0])

    def browse_dataset_loc(self):
        fname = QFileDialog.getOpenFileName(self, 'Open file', '', 'Datasets binary (*.exe, *)')
        if fname[0]:
            self.dataset_loc_input.setText(fname[0])

    def run_script(self):
        primers = self.primer_input.text().strip()
        primer_table = self.primer_table_input.text().strip()
        prefix = self.prefix_input.text().strip()
        working_directory = self.working_dir_input.text().strip()
        sassy_loc = self.sassy_input.text().strip()
        ncbi_virus = self.ncbi_virus_input.text().strip()
        ncbi_bacteria = self.ncbi_bacteria_input.text().strip()
        year = self.year_input.text().strip()
        directory_db = self.directory_db_input.text().strip()
        fasta_db = self.fasta_db_input.text().strip()
        glob_db = self.glob_db_input.text().strip()
        max_indel = int(self.max_indel_input.text().strip())
        max_mismatch = int(self.max_mismatch_input.text().strip())
        indel_mult = int(self.indel_mult_input.text().strip())
        threads = int(self.threads_input.text().strip())
        max_primer_dist = int(self.max_primer_dist_input.text().strip())
        dataset_loc = self.dataset_loc_input.text().strip()
        only_product_forming = self.only_product_forming_checkbox.isChecked()

        if not (primers or primer_table):
            QMessageBox.warning(self, 'Input Error', 'Please provide either a primer file or primer table.')
            return
        if primers and primer_table:
            QMessageBox.warning(self, 'Input Error', 'Please provide either a primer file OR primer table, not both.')
            return
        if not prefix or not working_directory:
            QMessageBox.warning(self, 'Input Error', 'Please provide prefix and working directory.')
            return
        if not ncbi_virus and not ncbi_bacteria and not directory_db and not fasta_db and not glob_db:
            QMessageBox.warning(self, 'Input Error', 'Please provide at least one genome database.')
            return

        # Get primer information
        if primers:
            primer_dict = get_primer_sequences(primers)
            primer_info = {}
        else:
            primer_dict, primer_info = get_primer_table(primer_table)

        if not os.path.exists(working_directory):
            os.makedirs(working_directory)
        elif os.path.exists(working_directory) and not os.path.isdir(working_directory):
            QMessageBox.critical(self, 'Error', 'An error occurred: working directory is not a directory.')
            return

        # Build database from all sources
        from collections import defaultdict
        amb_bases = defaultdict(set)
        
        primer_file = os.path.join(working_directory, prefix + "_db.fasta")
        open(primer_file, 'w').close()
        
        if ncbi_virus:
            fasta_files = download_virus(ncbi_virus, working_directory, prefix, date=year, datasets=dataset_loc)
            for fasta_file in fasta_files:
                amb_bases = get_db_fasta(fasta_file, working_directory, prefix, amb_bases)
        if ncbi_bacteria:
            fasta_files = download_bac(ncbi_bacteria, working_directory, prefix, datasets=dataset_loc)
            amb_bases = get_db_fastas(fasta_files, working_directory, prefix, amb_bases)
        if directory_db:
            fasta_files = get_db_folder(directory_db)
            amb_bases = get_db_fastas(fasta_files, working_directory, prefix, amb_bases)
        if fasta_db:
            amb_bases = get_db_fasta(fasta_db, working_directory, prefix, amb_bases)
        if glob_db:
            fasta_files = get_db_glob(glob_db)
            amb_bases = get_db_fastas(fasta_files, working_directory, prefix, amb_bases)

        self.worker = BlastWorker(primer_file, primer_dict, primer_info, working_directory, prefix, sassy_loc, max_indel, max_mismatch, indel_mult, threads, max_primer_dist)
        self.worker.finished.connect(self.on_blast_finished)
        self.worker.start()

        self.animation_label = QLabel('Running...')
        self.animation_label.setStyleSheet("background-color: lightgreen; font-weight: bold")
        self.layout().addWidget(self.animation_label, 13, 0, 1, 4)
        self.timer = QTimer(self)
        self.timer.timeout.connect(self.update_animation)
        self.timer.start(500)

    def update_animation(self):
        current_text = self.animation_label.text()
        if current_text.endswith('...'):
            self.animation_label.setText('Running')
        else:
            self.animation_label.setText(current_text + '.')

    def on_blast_finished(self):
        self.timer.stop()
        self.animation_label.setText('Finished')
        QMessageBox.information(self, 'Success', 'Script ran successfully.')

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = PrimerCheqGUI()
    sys.exit(app.exec_())   