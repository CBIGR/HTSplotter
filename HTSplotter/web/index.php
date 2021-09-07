<!doctype html>
<html lang="en">
	<head>
	<!-- Required meta tags -->
		<meta charset="utf-8">
		<meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">

		<!-- Bootstrap CSS -->
		<link rel="stylesheet" href="bootstrap-4.5.2-dist/css/bootstrap.min.css">
		<!-- own CSS -->
		<link rel="stylesheet" href="css/style.css">

		<title>HTSplotter - Home</title>
	</head>
	<body class="d-flex flex-column h-100">
		<header>
				<nav class="navbar navbar-expand-md navbar-dark bg-dark fixed-top">
				<a class="navbar-brand" href="index.php">HTSplotter</a>
				<button class="navbar-toggler" type="button" data-toggle="collapse" data-target="#navbarsExampleDefault" aria-controls="navbarsExampleDefault" aria-expanded="false" aria-label="Toggle navigation">
					<span class="navbar-toggler-icon"></span>
				</button>

				<div class="collapse navbar-collapse justify-content-end" id="navbarsExampleDefault">
					<ul class="nav">
						<li class="nav-item">
							<a class="nav-link" href="analysis.php">Run analysis</a>
						</li>
						<li class="nav-item">
							<a class="nav-link" href="images/HTSplotterManual.pdf" target="_blank">Manual</a>
						</li>
						<li class="nav-item">
							<a class="nav-link" href="contact.php">Feedback</a>
						</li>
					</ul>
				</div>
			</nav>
		</header>

		<main role="main" class="container flex-shrink-0">
			<div class="container-fluid">
                <h1 class="center">HTSplotter</h1>
                <h5>An end-to-end data processing, analysis and visualisation tool for chemical and genetic <i> in vitro </i>
                    perturbation screens</h5>
                <div class="row">&nbsp;</div>
				<div class="row">
					<div class="col-12">
						<div class="card">
							<div class="card-body">
								<h5 class="card-title">Example files</h5>
								<p class="card-text">Input example files for each experiment type. For more information, please check the <a href="images/HTSplotterManual.pdf" download rel="noopener noreferrer" target="_blank"> manual.</a></p>
								<img src="images/Mainplots_website2.png" alt="" border=3 width=100%></img>
								<table class="table table-responsive">
									<thead>
									<tr>
										<th scope="col">Experiment type</th>
										<th scope="col">Input file</th>
										<th scope="col">Results visualization</th>
										<th scope="col">Output files</th>
									</tr>
									</thead>
									<tbody>
									<tr>
										<td rowspan="3">Drug screen</td>
										<td><a href="demo/drugscreen_1timepoint.txt" download rel="noopener noreferrer" scope="col"> drugscreen_1timepoint</td>
										<td><a href="demo/results_pdf/drugscreen_1timepoint.pdf" download rel="noopener noreferrer" scope="col">pdf</td>
										<td><a href="demo/results_zip/drugscreen_1timepoint.zip" download rel="noopener noreferrer" scope="col">download</td>
									</tr>
									<tr>
										<td><a href="demo/drugscreen_severaltimepoint_1control.txt" download rel="noopener noreferrer" class="not-first-cell">drugscreen_severaltimepoint_1control</td>
										<td><a href="demo/results_pdf/drugscreen_severaltimepoint_1control.pdf" download rel="noopener noreferrer" scope="col">pdf</td>
										<td><a href="demo/results_zip/drugscreen_severaltimepoint_1control.zip" download rel="noopener noreferrer" scope="col">download</td>
									</tr>
									<tr>
										<td><a href="demo/drugscreen_severaltimepoint_severalcontrol.txt" download rel="noopener noreferrer" class="not-first-cell">drugscreen_severaltimepoint_severalcontrol</td>
										<td><a href="demo/results_pdf/drugscreen_severaltimepoint_severalcontrol.pdf" download rel="noopener noreferrer" scope="col">pdf</td>
										<td><a href="demo/results_zip/drugscreen_severaltimepoint_severalcontrol.zip" download rel="noopener noreferrer" scope="col">download</td>
									</tr>
									<tr>
										<td rowspan="3">Drug combination</td>
										<td><a href="demo/drug_combination_screen_1timepoint.txt" download rel="noopener noreferrer" class="not-first-cell">drug_combination_screen_1timepoint</td>
										<td><a href="demo/results_pdf/drug_combination_screen_1timepoint.pdf" download rel="noopener noreferrer" scope="col">pdf</td>
										<td><a href="demo/results_zip/drug_combination_screen_1timepoint.zip" download rel="noopener noreferrer" scope="col">download</td>
									</tr>
									<tr>
										<td><a href="demo/drug_combination_several_time_points.txt" download rel="noopener noreferrer" class="not-first-cell">drug_combination_several_time_points</td>
										<td><a href="demo/results_pdf/drug_combination_several_time_points.pdf" download rel="noopener noreferrer" scope="col">pdf</td>
										<td><a href="demo/results_zip/drug_combination_several_time_points.zip" download rel="noopener noreferrer" scope="col">download</td>
									</tr>
									<tr>
										<td><a href="demo/drug_combination_several_time_points_repetitive_conditions.txt" download rel="noopener noreferrer" class="not-first-cell">drug_combination_several_time_points_repetitive_conditions</td>
										<td><a href="demo/results_pdf/drug_combination_several_time_points_repetitive_conditions.pdf" download rel="noopener noreferrer" scope="col">pdf</td>
										<td><a href="demo/results_zip/drug_combination_several_time_points_repetitive_conditions.zip" download rel="noopener noreferrer" scope="col">download</td>
									</tr>
									<tr>
										<td rowspan="2">Genetic perturbagen</td>
										<td><a href="demo/gene_perturbagen_1timepoint_1control.txt" download rel="noopener noreferrer" class="not-first-cell">gene_perturbagen_1timepoint_1control</td>
										<td><a href="demo/results_pdf/gene_perturbagen_1timepoint_1control.pdf" download rel="noopener noreferrer" scope="col">pdf</td>
										<td><a href="demo/results_zip/gene_perturbagen_1timepoint_1control.zip" download rel="noopener noreferrer" scope="col">download</td>
									</tr>
									<tr>
										<td><a href="demo/gene_perturbagen_severaltimepoints.txt" download rel="noopener noreferrer" class="not-first-cell">gene_perturbagen_severaltimepoints</td>
										<td><a href="demo/results_pdf/gene_perturbagen_severaltimepoints.pdf" download rel="noopener noreferrer" scope="col">pdf</td>
										<td><a href="demo/results_zip/gene_perturbagen_severaltimepoints.zip" download rel="noopener noreferrer" scope="col">download</td>
									</tr>
									<tr>
										<td rowspan="2">Genetic-chemical perturbagen</td>
										<td><a href="demo/genetic-chemical_perturbagen_1time_point.txt" download rel="noopener noreferrer" class="not-first-cell">genetic-chemical_perturbagen_1time_point</td>
										<td><a href="demo/results_pdf/genetic-chemical_perturbagen_1time_point.pdf" download rel="noopener noreferrer" scope="col">pdf</td>
										<td><a href="demo/results_zip/genetic-chemical_perturbagen_1time_point.zip" download rel="noopener noreferrer" scope="col">download</td>
									</tr>
									<tr>
										<td><a href="demo/genetic-chemical_perturbagen_several-time_points.txt" download rel="noopener noreferrer" class="not-first-cell">genetic-chemical_perturbagen_several-time_points</td>
										<td><a href="demo/results_pdf/genetic-chemical_perturbagen_several-time_points.pdf" download rel="noopener noreferrer" scope="col">pdf</td>
										<td><a href="demo/results_zip/genetic-chemical_perturbagen_several-time_points.zip" download rel="noopener noreferrer" scope="col">download</td>
									</tr>
									</tr>
									</tbody>
								</table>
							</div>
						</div>
					</div>
				</div>
                <div class="row">&nbsp;</div>
                <div class="row">
                    <div class="col-12">
                        <h5 class="card-title">How to refer to this software?</h5>
                        <p class="card-text">Carolina Nunes, Jasper Anckaert, Fanny De Vloed, Jolien De Wyn,
                            Kaat Durinck, Jo Vandesompele, Frank Speleman and Vanessa Vermeirssen,
                            "Automatic end-to-end analysis of high-throughput<i> in vitro </i>
                            cell culture screening by HTSplotter", 2021 </p>
                        <p class="card-text">Source code available at
                            <a href="https://github.ugent.be/vermeirssenlab/HTSplotter">HTSplotter</a>, under GPLV3 license</p>
                    </div>
                </div>
                <div class="row">&nbsp;</div>
                <div class="row">&nbsp;</div>
                <div class="row">&nbsp;</div>
                <div class="row">&nbsp;</div>
                <div class="row">&nbsp;</div>
				<div class="row">&nbsp;</div>
			</div>
		</main>

		<footer class="footer mt-auto py-3 flex-shrink-0">	
			 <div class="container">
			 	<div class="row">
			 		<div class="col align-self-center">
			 			<a href="https://www.irc.ugent.be/index.php?id=vermeirssenhome" target="_blank" title="CBIGR"><img src="images/CBIGRlogo2.png" class="rounded footer_img" alt="CBIGR"></a>
						<a href="https://www.spelemanlab.org/" target="_blank" title="SpelemanLab"><img src="images/FSPlogo.png" class="rounded footer_img" alt="SpelemanLab"></a>
						<a href="https://oncornalab.ugent.be" target="_blank" title="OncoRNAlab"><img src="images/oncolab_1.png" class="rounded footer_img" alt="OncoRNAlab"></a>
						<a href="https://www.ugent.be/en" target="_blank" title="UGent"><img src="images/UGent.png" class="rounded footer_img" alt="UGent"></a>
			 		</div>
				</div>
			</div>
		</footer>

		<!-- Optional JavaScript -->
		<!-- jQuery first, then Popper.js, then Bootstrap JS -->
		<script src="https://code.jquery.com/jquery-3.6.0.min.js" integrity="sha256-/xUj+3OJU5yExlq6GSYGSHk7tPXikynS7ogEvDej/m4=" crossorigin="anonymous"></script>
		<script src="https://cdn.jsdelivr.net/npm/popper.js@1.16.1/dist/umd/popper.min.js" integrity="sha384-9/reFTGAW83EW2RDu2S0VKaIzap3H66lZH81PoYlFhbGU+6BZp6G7niu735Sk7lN" crossorigin="anonymous"></script>
		<script src="bootstrap-4.5.2-dist/js/bootstrap.min.js"></script>
		<script src="js/action.js"></script>
	</body>
</html>
