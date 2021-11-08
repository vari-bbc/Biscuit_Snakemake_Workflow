# Creating Figures of Included DAGs
For those who want to see a figure of the workflows for either the one or five sample example workflows, the following
commands will show you how to create these:
```
# One sample files
dot -Tpng one_sample_DAG_default_workflow > one_sample_DAG_default_workflow.png
dot -Tpng one_sample_DAG_work_workflow > one_sample_DAG_work_workflow.png

# Five sample files
dot -Tpng five_sample_DAG_default_workflow > five_sample_DAG_default_workflow.png
dot -Tpng five_sample_DAG_work_workflow > five_sample_DAG_work_workflow.png
```

The default workflows were generated using the default `config.yaml` file, while the full workflows set all optional
rules to `True` before creating the DAG file (`snakemake --dag > DAG_workflow`).
