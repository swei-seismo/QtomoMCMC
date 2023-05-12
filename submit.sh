output=$(sbatch run.sb)
echo "$output"
prev_job_id=$(echo "$output" | grep -oP 'Submitted batch job \K\d+')
sbatch --qos=scavenger --dependency=afterok:$prev_job_id Plot_in_GMT.sb

