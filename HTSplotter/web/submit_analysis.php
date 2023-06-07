<?php
include_once('functions.php');

if ($_POST['user_input'] == 0){
	$date_part = date("Ymd");
	$user_string = generateRandomString();
	$cmd_dir = 'user_files/User_'.$user_string.'_'.$date_part.'_HTSplotter/';
	$user_dir = $cmd_dir.'Inputfile/';
	$response_array['user_input'] = $cmd_dir;

	$oldumask = umask(0); 
	mkdir($user_dir, 0777, TRUE);
	umask($oldumask);

	$total = count($_FILES['csv_file']['name']);

	for ($i=0; $i < $total; $i++ ) {
		$target_file = $user_dir . basename($_FILES['csv_file']['name'][$i]);
		move_uploaded_file($_FILES['csv_file']['tmp_name'][$i], $target_file);
	}


	$command = 'python3 /var/www/html/scripts/web.py -f /var/www/html/'.$cmd_dir.' -b '.escapeshellarg($_POST['bio_rep']).' -n '.escapeshellarg($_POST['biorep_fn']).' -e '.escapeshellarg($_POST['exp_effect']).' -i '.escapeshellarg($_POST['info_readout']).' -r '.escapeshellarg($_POST['readout_unit']).' -sy '.escapeshellarg($_POST['syn_ana']).' -u '.escapeshellarg('no').' -sh '.escapeshellarg($_ENV['S3_HOST']).' -sb '.escapeshellarg($_ENV['S3_BUCKET']).' -sa '.escapeshellarg($_ENV['S3_ACCESS_KEY']).' -ss '.escapeshellarg($_ENV['S3_SECRET_KEY']);
	$output = shell_exec($command);

	$errorfile = '';
	$error_array = glob($user_dir.'*_Errorfile.txt');
	if (!empty($error_array)) {
		$errorfile = $error_array[0];
	}
	if (file_exists($errorfile) && filesize($errorfile) > 0) {
		$response_array["errorfile"] = $errorfile;
	} else {
		$oparray = preg_split("#[\r\n]+#", trim($output));
		$exp_types = preg_grep("/^Experiment type/i", $oparray);
		$exp_types = array_values($exp_types);
		$filenames = preg_grep("/^File/i", $oparray);
		$filenames = array_values($filenames);
		$exptype_str = '';
		$filename_str = '';
		$name = array_values($filenames);
		$name_str = '';
		$infofile = '';
		foreach ($exp_types as $key => $exp_type) {
			$exptype_str = $exptype_str.trim(explode(": ", $exp_types[$key])[1]).",";
			$filename_str = $filename_str.trim(explode(": ", $filenames[$key])[1]).",";
			$name_str1 = $name_str.trim(explode(": ", $name[$key])[1]).",";
			$name_str2 = explode(".", $name_str1);
			$infofile = $infofile.trim($cmd_dir.$name_str2[0].'_information.txt').",";
		}
		$response_array["exp_types"] = substr($exptype_str, 0, -1);
		$response_array["filenames"] = substr($filename_str, 0, -1);
		$response_array["infofile"] = substr($infofile, 0, -1);
	}
	echo json_encode($response_array);
} else {
	$exptype_str = "'".str_replace(",","' '",$_POST['exp_types'])."'";
	$command = 'python3 /var/www/html/scripts/web.py -f /var/www/html/'.$_POST['user_input'].' -b '.escapeshellarg($_POST['bio_rep']).' -n '.escapeshellarg($_POST['biorep_fn']).' -e '.escapeshellarg($_POST['exp_effect']).' -i '.escapeshellarg($_POST['info_readout']).' -r '.escapeshellarg($_POST['readout_unit']).' -sy '.escapeshellarg($_POST['syn_ana']).' -u '.escapeshellarg('yes').' -t '.$exptype_str.' -sh '.escapeshellarg($_ENV['S3_HOST']).' -sb '.escapeshellarg($_ENV['S3_BUCKET']).' -sa '.escapeshellarg($_ENV['S3_ACCESS_KEY']).' -ss '.escapeshellarg($_ENV['S3_SECRET_KEY']);
	$output = shell_exec($command);
	$oparray = preg_split("#[\r\n]+#", trim($output));
	$result = preg_grep("/^Output/i", $oparray);
	$result = array_values($result);
	$response_array['pdf'] = trim(explode(": ", $result[0])[1]);
	$response_array['html'] = str_replace('/var/www/html/', '', glob('/var/www/html/'.$_POST['user_input'].'Output_html/' . "{*.html}", GLOB_BRACE));
	$response_array['user_input'] = $_POST['user_input'];
	$response_array['zip'] = 'https://'.$_ENV['S3_HOST'].'/'.$_ENV['S3_BUCKET'].'/'.explode("/", $_POST['user_input'])[1].'.zip';
	//$response_array['zip'] = 'https://'.$_ENV['S3_HOST'].'/'.$_ENV['S3_BUCKET'].'/'.explode("/", $_POST['user_input'])[1].'.zip';

	echo json_encode($response_array);
}

?>
