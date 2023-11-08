// This Jenkinsfile is used by Jenkins to run the 'UniProt Update' step of Reactome's release.
// This step synchronizes Reactome's database with UniProt.

import org.reactome.release.jenkins.utilities.Utilities

// Shared library maintained at 'release-jenkins-utils' repository.
def utils = new Utilities()

pipeline {
	agent any

	stages {
		// This stage checks that an upstream step, ConfirmReleaseConfigs, was run successfully.
		stage('Check ConfirmReleaseConfigs build succeeded'){
			steps{
				script{
					utils.checkUpstreamBuildsSucceeded("ConfirmReleaseConfigs")
				}
			}
		}
		// Download uniprot_sprot.xml.gz from UniProt.
		stage('Setup: Download uniprot_sprot.xml.gz'){
			steps{
				script{
					sh "wget -N ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.xml.gz"
				}
			}
		}
		// This stage backs up the gk_central database before it is modified.
		stage('Setup: Back up gk_central before modifications'){
			steps{
				script{
					withCredentials([usernamePassword(credentialsId: 'mySQLCuratorUsernamePassword', passwordVariable: 'pass', usernameVariable: 'user')]){
					    utils.takeDatabaseDumpAndGzip("${env.GK_CENTRAL_DB}", "uniprot_update", "before", "${env.CURATOR_SERVER}")
					}
				}
			}
		}
		// This stage executes the UniProt Update Java program.
		stage('Main: UniProt Update'){
			steps {
				script{
					withCredentials([file(credentialsId: 'Config', variable: 'ConfigFile')]){
						sh "java -Xmx${env.JAVA_MEM_MAX}m -jar target/uniprot*-jar-with-dependencies.jar $ConfigFile"
					}
				}
			}
		}
		// This stage backs up the gk_central database after modification.
		stage('Post: Backup gk_central after modifications'){
			steps{
				script{
					withCredentials([usernamePassword(credentialsId: 'mySQLCuratorUsernamePassword', passwordVariable: 'pass', usernameVariable: 'user')]){
						utils.takeDatabaseDumpAndGzip("${env.GK_CENTRAL_DB}", "uniprot_update", "after", "${env.CURATOR_SERVER}")
					}
				}
			}
		}
		// This stage emails the uniprot.wiki file to the default recipients list.
		stage('Post: Email UniProt.wiki file'){
			steps{
				script{
					def uniprotWikiFile = "uniprot.wiki"

					def releaseVersion = utils.getReleaseVersion();
					def emailSubject = "UniProt Update Reports for v${releaseVersion}"
					def emailBody = "Hello,\n\nThis is an automated message from Jenkins regarding an update for v${releaseVersion}. The UniProt Update step has completed. Please review the ${uniprotWikiFile} file attached to this email. If it looks correct, the contents of the file need to be uploaded to https://devwiki.reactome.org/index.php/Reports_Archive under 'UniProt Update Reports'. Please add the current UniProt wiki URL to the 'Archived reports' section of the page. If the file looks incorrect, please email the developer running Release. \n\nThanks!"
					utils.sendEmailWithAttachment("$emailSubject", "$emailBody", "$uniprotWikiFile")
				}
			}
		}
		// All databases, logs, and data files generated by this step are compressed before moving them to the Reactome S3 bucket.
		// All files are then deleted. TODO: Once this step has changed to a Java module, the archive module from the shared library should be used.
		stage('Post: Archive Outputs'){
			steps{
				script{
				    def releaseVersion = utils.getReleaseVersion();
					def s3Path = "${env.S3_RELEASE_DIRECTORY_URL}/${releaseVersion}/uniprot_update"
					sh "mkdir -p databases/ data/ logs/"

					sh "mv --backup=numbered *_${releaseVersion}_*.dump.gz databases/"
					sh "mv ewasCoordinatesReport.txt data/"
					sh "mv sequence_uniprot_report.txt data/"
					sh "mv reference_DNA_sequence_report.txt data/"
					sh "mv duplicated_db_id.txt data/"
					sh "mv trembl_to_update.acc data/"
					sh "mv uniprot_sprot.xml data/"
					sh "mv uniprot.wiki data/"
					sh "mv uniprot.* logs/"

					sh "gzip -r data/* logs/*"
					sh "aws s3 --no-progress --recursive cp databases/ $s3Path/databases/"
					sh "aws s3 --no-progress --recursive cp logs/ $s3Path/logs/"
					sh "aws s3 --no-progress --recursive cp data/ $s3Path/data/"
					sh "rm -r databases logs data"
				}
			}
		}
	}
}