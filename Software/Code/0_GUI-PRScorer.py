from bottle import Bottle, request, run, static_file, template, response
import os, time, subprocess

app = Bottle()

# Chemins relatifs à ce fichier (0_GUI-PRScorer.py)
#Relative paths to this file (0_GUI-PRScorer.py)
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
PIPELINE_DIR = os.path.abspath(os.path.join(BASE_DIR, "../../"))  # Va de Code/ → Pipeline/ #Goes from Code/ → Pipeline/

DATA_DIR = os.path.join(PIPELINE_DIR, "Data")
RESULTS_DIR = os.path.join(PIPELINE_DIR, "Results")
VIEWS_DIR = os.path.join(BASE_DIR, "Views")  # index.tpl (code HTML)
CODE_DIR = BASE_DIR  # shell scripts sont dans Code/
PROGRESS_FILE = os.path.join(BASE_DIR, "Progress.txt")  # Forgotten code to show progress,bottle does not support dynamic updates, so this file is not used in the current implementation

@app.route('/')
def home():
    return template("index", template_lookup=[VIEWS_DIR])

@app.post('/upload')
def upload_folder():
    folder_name = request.forms.get('folder_name')
    file_type_choice = request.forms.get('FileTypeChoice')
    bead_chip_choice = request.forms.get('BeadChipChoice')

    if not folder_name:
        return "Folder name is required."

    target_dir = os.path.join(DATA_DIR, folder_name, file_type_choice)
    os.makedirs(target_dir, exist_ok=True)

    upload = request.files.getall('files')
    for file in upload:
        save_path = os.path.join(target_dir, file.filename)
        file.save(save_path, overwrite=True)

    for file in upload:
        save_path = os.path.join(target_dir, file.filename)
        while True:
            initial_size = os.path.getsize(save_path)
            time.sleep(1)
            current_size = os.path.getsize(save_path)
            if initial_size == current_size:
                break

    if file_type_choice == "IDAT":
        subprocess.run(["sh", os.path.join(CODE_DIR, "1_idat2gtc.sh"), folder_name, bead_chip_choice])
    elif file_type_choice == "GTC":
        subprocess.run(["sh", os.path.join(CODE_DIR, "2_gtc2vcf.sh"), folder_name, bead_chip_choice])

    return f"Folder uploaded successfully to {target_dir} processing {file_type_choice} files."

@app.post('/Post_Imputation')
def post_imputation():
    folder_name = request.forms.get('folder_name')
    DL_link = request.forms.get('DL_Link')

    if not DL_link:
        return "Download link is required."

    subprocess.run(["sh", os.path.join(CODE_DIR, "5_Imputation_Download.sh"), DL_link, folder_name])
    return f'''The link "{DL_link}" was successfully posted.
        <br>
        <form action="/" method="get">
        <input type="submit" value="Run another analysis">
        </form>
    '''

@app.route('/progress')
def get_progress():
    if not os.path.exists(PROGRESS_FILE):
        return {'output': '[Progress file not found]'}
    with open(PROGRESS_FILE, 'r', encoding='utf-8') as f:
        return {'output': f.read()}

@app.route('/results_list')
def get_results_list():
    sample_list_path = os.path.join(RESULTS_DIR, 'SampleList.txt')
    if not os.path.exists(sample_list_path):
        response.status = 404
        return "SampleList.txt not found"
    with open(sample_list_path, 'r', encoding='utf-8') as f:
        return {'samples': [line.strip() for line in f if line.strip()]}

@app.route('/results/<item>')
def get_results(item):
    folder = os.path.join(RESULTS_DIR, item)
    if not os.path.exists(folder):
        response.status = 404
        return "No results found for this sample."

    result_texts, images = [], []

    for filename in os.listdir(folder):
        full_path = os.path.join(folder, filename)
        rel_path = f"/results_file/{item}/{filename}"
        if filename.endswith('.txt'):
            with open(full_path, 'r', encoding='utf-8') as f:
                result_texts.append(f.read())
        elif filename.endswith(('.svg', '.png')):
            images.append(rel_path)

    # Trier les images : celles finissant par PC1_PC2.svg d'abord puis celles finissant par PC2_PC3.svg
    # Sort images: those ending with PC1_PC2.svg first then PC2_PC3.svg
    def sort_key(path):
        if path.endswith("PC1_PC2.svg"):
            return (0, path)
        elif path.endswith("PC2_PC3.svg"):
            return (1, path)
        else:
            return (2, path)

    images.sort(key=sort_key)

    return {'texts': result_texts, 'images': images}


@app.route('/results_file/<item>/<filename:path>')
def serve_result_file(item, filename):
    return static_file(filename, root=os.path.join(RESULTS_DIR, item))

if __name__ == '__main__':
    run(app, host='0.0.0.0', port=8080, debug=True)
