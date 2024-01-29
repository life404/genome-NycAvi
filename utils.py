#! /usr/bin/python3

import subprocess
import concurrent.futures as cf
from pathlib import Path
from tqdm import tqdm
import argparse
import shutil


def cmd_run(command, ou_path):
    run = subprocess.Popen(
        command,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        cwd=ou_path,
    )
    output = run.communicate()
    return_code = run.wait()
    return return_code, output, command


def cmd_run_multiple(command_lst, ou_path, check, verbose=False, num=5, quite=False):
    check_file = Path(ou_path) / f"{check}.ok"
    error_log = Path(ou_path) / f"{check}.error"
    if error_log.exists():
        error_log.unlink()
    if isinstance(command_lst, str):
        command_lst = [command_lst]

    if check_file.exists():
        print(f"FIND THE CHECK POINT `{check_file.name}`, SKIP THIS STEP")
    else:
        print(f"RUN THOES COMMAND IN PARALLEL:")
        print("+" * 64)
        if len(command_lst) > 5:
            for command in command_lst[1:6]:
                print(command)
            print(
                "\ntoo much commands, the commands are more than 10, only show the top 5 commands ......\n"
            )
        else:
            for command in command_lst:
                print(command)

        num = min(len(command_lst), num)
        print("The Progressing:\n")
        progress_bar = tqdm(total=len(command_lst))
        with cf.ProcessPoolExecutor(max_workers=num) as e:
            process_lst = [
                e.submit(cmd_run, command, ou_path) for command in command_lst
            ]
            return_code_lst = list()
            output_lst = list()
            error_lst = list()
            for process in cf.as_completed(process_lst):
                return_code, output, command = process.result()
                return_code_lst.append(return_code)
                if return_code != 0:
                    error_lst.append(output[0])
                    error_lst.append(command + "\n")
                output_lst.append(output)
                progress_bar.update(1)

        progress_bar.close()

        if verbose:
            for output in output_lst:
                for line in output:
                    print(line)

        if sum(return_code_lst) == 0:
            check_file.touch()
            if check == 'null':
                check_file.unlink()
            print("+" * 64)
            print(f"JOB(s) has FINISHED: the check file is {check}")
        else:
            print("+" * 64)
            print(
                f"SOMETHING WRONG, PLEASE SET `verbose = True` or SEE {error_log} TO CHECK THE WRROR"
            )
            # print(error_lst)
            error_content = open(error_log, "w")
            error_content.writelines(error_lst)
            error_content.close()
            exit(1)
    return check_file


def cmd_run_R(rconsole, rcommand, ou_path, check):
    check = Path(ou_path) / f"{check}.ok"
    if not isinstance(rcommand, str):
        print("The R command run in R cosole should be a str")
        exit(1)

    if check.exists():
        print(f"FIND THE CHECK POINT `{check.name}`, SKIP THIS STEP")
    else:
        print("RUN THOSE R COMMANDS:\n")
        print("+" * 64)
        print(rcommand)
        rconsole(rcommand)
        check.touch()
        print("+" * 64)
        print(f"JOB(s) has FINISHED: the check file is {check}")

    return check
