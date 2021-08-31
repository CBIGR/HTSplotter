<?php

include 'config.php';

$post = (!empty($_POST)) ? true : false;

if($post)
{

	$name = stripslashes($_POST['fname']).' '.stripslashes($_POST['lname']);
	$email = trim($_POST['email']);
	$subject = "[HTSplotter] - Message/question";
	$message = stripslashes($_POST['message']);


	$error = '';



	if(!$error)
	{
		$mail = mail(WEBMASTER_EMAIL, $subject, $message,
		     "From: ".$name." <".$email.">\r\n"
		    ."Reply-To: ".$email."\r\n"
		    ."Bcc: ".BCC_EMAIL."\r\n"
		    ."X-Mailer: PHP/" . phpversion());


		if($mail)
		{
			echo 'OK';
		}

	}


}
else {
	echo 'NOT OK';
}
?>