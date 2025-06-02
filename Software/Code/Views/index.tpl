<!DOCTYPE html>
<html lang="fr">
<head>
    <meta charset="UTF-8">
    <title>PRScorer</title>
        <style>
        body {
            margin: 0;
            font-family: Arial, sans-serif;
            background-color: #f7f9fc;
            color: #333333;
            display: flex;
            flex-direction: column;
            align-items: center;
            padding: 2em;
        }

        h1, h2 {
            text-align: center;
            color: #333;
        }

        form {
            background-color: #ffffff;
            border-radius: 10px;
            box-shadow: 0 2px 10px rgba(0, 0, 0, 0.05);
            padding: 2em;
            width: 100%;
            max-width: 600px;
            margin-bottom: 2em;
        }

        label {
            display: block;
            margin-top: 1em;
            font-weight: bold;
        }

        input[type="text"],
        select,
        input[type="file"] {
            width: 100%;
            padding: 0.75em;
            margin-top: 0.5em;
            border: 1px solid #ccc;
            border-radius: 5px;
            background-color: #f9f9f9;
        }

        input[type="submit"] {
            margin-top: 1.5em;
            background-color: #007bff;
            color: #ffffff;
            border: none;
            padding: 0.75em 1.5em;
            font-size: 1em;
            border-radius: 5px;
            cursor: pointer;
            transition: background-color 0.3s ease;
        }

        input[type="submit"]:hover {
            background-color: #0056b3;
        }

        .box {
            background-color: #ffffff;
            border: 1px solid #ccc;
            border-radius: 10px;
            padding: 1em;
            margin-top: 1em;
            white-space: pre-wrap;
            width: 100%;
            max-width: 600px;
        }

        .dropdown {
            position: relative;
            width: 100%;
            max-width: 600px;
        }

        .dropdown-content {
            display: none;
            position: absolute;
            background-color: #ffffff;
            border: 1px solid #ccc;
            border-top: none;
            z-index: 1000;
            max-height: 200px;
            overflow-y: auto;
            width: 100%;
        }

        .dropdown-content div {
            padding: 0.75em;
            cursor: pointer;
        }

        .dropdown-content div:hover {
            background-color: #f0f0f0;
        }

        #searchInput {
            width: 100%;
            padding: 0.75em;
            margin-top: 0.5em;
            border: 1px solid #ccc;
            border-radius: 5px;
        }

        #resultsImages img {
            max-width: 100%;
            display: block;
            margin: 1em 0;
        }

        #folder_name_error {
            color: red;
            display: none;
            font-size: 0.9em;
        }
        .spinner {
        border: 2px solid #f3f3f3;
        border-top: 2px solid #ffffff;
        border-radius: 50%;
        width: 16px;
        height: 16px;
        animation: spin 0.8s linear infinite;
        display: inline-block;
        vertical-align: middle;
        margin-right: 8px;
        }

        @keyframes spin {
            0% { transform: rotate(0deg); }
            100% { transform: rotate(360deg); }
        }
    </style>
</head>
<body>
    <h1>Welcome to the PRScorer</h1>

    <!-- Raw file transformation and imputation part -->
    <h2>Upload some files</h2>
    <form action="/upload" method="post" enctype="multipart/form-data" onsubmit="return validateFolderName()">
        <label>Analysis Name:</label>
        <input type="text" id="folder_name" name="folder_name" required oninput="checkFolderName()">
        <span id="folder_name_error" style="color: red; display: none;">No space or special characters allowed.</span>
        <br><br>
        <label>Data type:</label>
        <select id="FileType" name="FileTypeChoice" onchange="updateAccept()">
            <option value="IDAT">idat</option>
            <option value="GTC" selected>gtc</option>
        </select>
        <br><br>
        <label>BeadChip type:</label>
        <select id="BeadChipType" name="BeadChipChoice" required>
            <option value="Infinium_Omni5-4_v1.2">Infinium Omni5-4v1.2</option> <!-- Values should be the EXACT name as the file containing the cluster and manifest -->
            <option value="Global_Diversity_Array_with_PRS">Global Diversity Array with PRS</option>
        </select>
        <br><br>
        <input type="file" id="files" name="files" accept=".gtc" multiple required>
        <br><br>
        <input type="submit" id="uploadBtn" value="Upload">
    </form>

    <!-- Boîte de progression -->
    <!-- <div class="box" id="progressBox">[Waiting for progress...]</div>-->

    <hr>

    <!-- Post-imputation part -->
    <h2>I have received an imputation email</h2>
    <form action="/Post_Imputation" method="post">
        <label>Analysis Name:</label>
        <input type="text" name="folder_name" required>
        <br><br>
        <label>CURL download link:</label>
        <input type="text" name="DL_Link" required>
        <br><br>
        <input type="submit" id="postBtn" value="Post_Imputation">
    </form>
    <hr>
    <!--  View Results avec barre de recherche -->
    <h2>View Results</h2>
    <div class="dropdown">
        <input type="text" id="searchInput" onkeyup="filterSamples()" placeholder="Search for a sample...">
        <div id="dropdownList" class="dropdown-content">
            <!-- Liste remplie dynamiquement -->
        </div>
    </div>

    <div class="box" id="resultsText">[Results will appear here]</div>
    <div id="resultsImages"></div>

    <!-- JavaScript -->
    <script>
        function updateAccept() {
            const type = document.getElementById('FileType').value;
            document.getElementById('files').accept = type === 'IDAT' ? '.idat' : '.gtc';
        }

        function checkFolderName() {
            const name = document.getElementById('folder_name').value;
            const regex = /^[a-zA-Z0-9_]+$/;
            document.getElementById('folder_name_error').style.display = regex.test(name) ? 'none' : 'inline';
        }

        function validateFolderName() {
            const regex = /^[a-zA-Z0-9_-]+$/;
            const name = document.getElementById('folder_name').value;
            return regex.test(name) || !alert("Invalid folder name");
        }

        //  Affichage dynamique de Progress.txt
        /*
        setInterval(() => {
            fetch('/progress')
                .then(res => res.json())
                .then(data => document.getElementById('progressBox').textContent = data.output);
        }, 5000);
        */

        // Dropdown + filtre
        let allSamples = [];

        fetch('/results_list')
            .then(res => res.json())
            .then(data => {
                allSamples = data.samples;
                updateDropdown(allSamples);
            });

        function updateDropdown(samples) {
            const list = document.getElementById("dropdownList");
            list.innerHTML = '';
            samples.forEach(sample => {
                const div = document.createElement("div");
                div.textContent = sample;
                div.onclick = () => {
                    document.getElementById("searchInput").value = sample;
                    list.style.display = "none";
                    loadResults(sample);
                };
                list.appendChild(div);
            });
        }

        function filterSamples() {
            const input = document.getElementById("searchInput").value.toLowerCase();
            const filtered = allSamples.filter(sample => sample.toLowerCase().includes(input));
            updateDropdown(filtered);
            document.getElementById("dropdownList").style.display = "block";
        }

        function loadResults(item) {
            fetch(`/results/${item}`)
                .then(res => res.json())
                .then(data => {
                    document.getElementById("resultsText").textContent = data.texts.join("\n\n");
                    const imgDiv = document.getElementById("resultsImages");
                    imgDiv.innerHTML = '';
                    data.images.forEach(src => {
                        const img = document.createElement('img');
                        img.src = src;
                        img.style.maxWidth = "600px";
                        img.style.display = "block";
                        img.style.marginBottom = "1em";
                        imgDiv.appendChild(img);
                    });
                });
        }

        // gestion ouverture/fermeture dropdown
        document.getElementById("searchInput").addEventListener("focus", () => {
            document.getElementById("dropdownList").style.display = "block";
        });

        document.addEventListener("click", function (e) {
            if (!e.target.closest('.dropdown')) {
                document.getElementById("dropdownList").style.display = "none";
            }
        });
        document.querySelector('form[action="/upload"]').addEventListener('submit', function () {
    const btn = document.getElementById("uploadBtn");
    btn.disabled = true;
    btn.innerHTML = `<span class="spinner"></span>Uploading...`;
});

document.querySelector('form[action="/Post_Imputation"]').addEventListener('submit', function () {
    const btn = document.getElementById("postBtn");
    btn.disabled = true;
    btn.innerHTML = `<span class="spinner"></span>Posting...`;
});

    </script>
</body>
</html>
