


#include <Python.h>
#include <silo.h>
#include <iostream>
#include <stdio.h>
#include <vector>


static PyObject * parse_file(PyObject *self, PyObject *args, PyObject *keywords)
{
  std::string ref_dir("domain_");

  //---------------------------------------------------------------------
  // Parse inputs
  //---------------------------------------------------------------------
  // Command line arguments
  const char *file_name = "../../../plot_687492";
  const char *field_type = "FaceFields";
  PyObject *field_names;
  const char *parallel_folder="";
  size_t npar=1;

  static char *keywordlist[] = {"file_name", "field_type", "field_names", "parallel_folder", "npar", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, keywords, "ssO|si", keywordlist, &file_name,
                                   &field_type, &field_names, &parallel_folder, &npar))
  {
    std::cout << "Error: inputs were not specified correctly!  Correct argments:" << std::endl;
    std::cout << "parse_file(file_name, field_type, field_names)" << std::endl;
    return NULL;
  }

  // Check the field_names input
  std::vector<std::string> field_names_requested;
  if (PyList_Check(field_names))
  {
    size_t field_names_size = PyList_Size(field_names);
    if (field_names_size > 0)
    {
      field_names_requested.resize(field_names_size);
      for (size_t ii=0; ii<field_names_size; ++ii)
      {
        PyObject *list_val = PyList_GetItem(field_names, ii);
        if (PyUnicode_Check(list_val))
        {
          field_names_requested[ii] = PyUnicode_AsUTF8(list_val);
        }
        else
        {
          std::cout << "Error: field_names[" << ii << "] is not a string!" << std::endl;
        }
      }
    }
    else
    {
      std::cout << "Error: field_names is empty" << std::endl;
    }
  }
  else
  {
    std::cout << "Error: field_names must be a list containing the field names to be exported!" << std::endl;
  }

  if (field_names_requested.size() == 0)
  {
    return NULL;
  }

  field_names_requested.push_back("FaceCenter");


  //---------------------------------------------------------------------
  // Check for available field values
  //---------------------------------------------------------------------
  // Open up the root silo
  std::vector<std::string> field_names_targets;
  std::vector<size_t> field_dimensions;
  int total_size = 0;

  DBfile *db = DBOpen(file_name, DB_UNKNOWN, DB_READ);
  DBtoc *toc = DBGetToc(db);
  size_t N = toc->ndir;

  for (size_t ii=0; ii<N; ++ii)
  {
    toc = DBGetToc(db);
    std::string dir(toc->dir_names[ii]);

    if (dir.find(ref_dir)!=std::string::npos)
    {
      // This is a domain directory
      DBSetDir(db, dir.c_str());
      toc = DBGetToc(db);
      size_t M = toc->ndir;

      for (size_t jj=0; jj<M; ++jj)
      {
        toc = DBGetToc(db);
        std::string sub_dir(toc->dir_names[jj]);

        if (sub_dir.find(field_type)!=std::string::npos)
        {
          // This is the field directory
          // Dive into it and test if any of the requested keys are available
          DBSetDir(db, sub_dir.c_str());
          toc = DBGetToc(db);

          for (int kk=0; kk<toc->nucdvar; ++kk)
          {
            for (size_t xx=0; xx<field_names_requested.size(); ++xx)
            {
              if (strcmp(toc->ucdvar_names[kk], field_names_requested[xx].c_str()) == 0)
              {
                field_names_targets.push_back(field_names_requested[xx]);
                DBucdvar *var = DBGetUcdvar(db, field_names_requested[xx].c_str());
                field_dimensions.push_back(var->nvals);
                total_size += var->nvals;
              }
            }
          }
          DBSetDir(db, "..");
        }
      }
      DBSetDir(db, "..");
      break;
    }
  }

  if (field_names_targets.size() == 0)
  {
    std::cout << "Error: no values to export were located" << std::endl;
    return NULL;
  }

  DBClose(db);


  //---------------------------------------------------------------------
  // Parse values
  //---------------------------------------------------------------------
  // Note: the last three positions are by default FaceCenter
  double current_time = 0.0;
  size_t current_cycle = 0;
  std::vector< std::vector<double> > export_values(total_size);

  // Process base and sub-files
  for (size_t ii=0; ii<npar; ++ii)
  {
    // Get current name
    std::string current_file = file_name;
    if (ii > 0)
    {
      int n;
      char buffer[100];
      n = sprintf(buffer, "%s/%s.%03d", parallel_folder, file_name, int(ii));
      current_file = buffer;
    }

    // Open the target file, and dive into the available domains
    DBfile *db = DBOpen(current_file.c_str(), DB_UNKNOWN, DB_READ);
    DBtoc *toc = DBGetToc(db);
    size_t N = toc->ndir;

    for (size_t jj=0; jj<N; ++jj)
    {
      // Save cycle, time
      if (ii == 0)
      {
        current_cycle = *((int *)DBGetVar(db, "cycle"));
        current_time = *((double *)DBGetVar(db, "dtime"));
      }

      toc = DBGetToc(db);
      std::string dir(toc->dir_names[jj]);

      if (dir.find(ref_dir)!=std::string::npos)
      {
        // Dive in and extract values
        DBSetDir(db, dir.c_str());
        DBSetDir(db, field_type);

        // Construct a filter (Only export values with ghostRank=0, flowFaceType=1)
        std::vector<int> export_indices;

        if (strcmp(field_type, "FaceFields") == 0)
        {
          std::vector<int> export_indices_tmp;
          DBucdvar *flowFaceType_var = DBGetUcdvar(db, "flowFaceType");
          void *flowFaceType_vals = flowFaceType_var->vals[0];

          for (int kk=0; kk<flowFaceType_var->nels; kk++)
          {
            int ff = ((int*)flowFaceType_vals)[kk];

            if (ff == 1)
            {
              export_indices_tmp.push_back(kk);
            }
          }

          DBucdvar *ghostRank_var = DBGetUcdvar(db, "ghostRank");
          void *ghostRank_vals = ghostRank_var->vals[0];

          for (size_t kk=0; kk<export_indices_tmp.size(); kk++)
          {
            int gr = ((int*)ghostRank_vals)[export_indices_tmp[kk]];
            if (gr < 0)
            {
              export_indices.push_back(export_indices_tmp[kk]);
            }
          }
        }
        else
        {
          DBucdvar *ghostRank_var = DBGetUcdvar(db, "ghostRank");
          void *ghostRank_vals = ghostRank_var->vals[0];

          for (int kk=0; kk<ghostRank_var->nels; kk++)
          {
            int gr = ((int*)ghostRank_vals)[kk];
            if (gr < 0)
            {
              export_indices.push_back(kk);
            }
          }
        }

        size_t new_values = export_indices.size();

        // Copy over data
        if (new_values > 0)
        {
          size_t xx = 0;
          for (size_t kk=0; kk<field_names_targets.size(); kk++)
          {
            DBucdvar *current_var = DBGetUcdvar(db, field_names_targets[kk].c_str());

            for (size_t yy=0; yy<field_dimensions[kk]; yy++)
            {
              size_t start_size = export_values[xx].size();
              export_values[xx].resize(start_size + new_values);
              void *current_vals = current_var->vals[yy];

              for (size_t zz=0; zz<new_values; zz++)
              {
                export_values[xx][zz + start_size] = ((double *)current_vals)[export_indices[zz]];
              }
              xx++;
            }
          }
        }

        DBSetDir(db, "../..");
      }
    }
    DBClose(db);
  }

  //---------------------------------------------------------------------
  // Build Outputs
  //---------------------------------------------------------------------
  // Setup the dictionary that will be returned
  PyObject *plot_dict = PyDict_New();
  PyDict_SetItem(plot_dict, PyUnicode_FromString("time"), PyFloat_FromDouble(current_time));
  PyDict_SetItem(plot_dict, PyUnicode_FromString("cycle"), PyLong_FromLong(current_cycle));

  size_t xx = 0;
  for (size_t ii=0; ii<field_names_targets.size(); ++ii)
  {
    for (size_t jj=0; jj<field_dimensions[ii]; ++jj)
    {
      std::string output_key = field_names_targets[ii];
      if (field_dimensions[ii] > 1)
      {
        int n;
        char buffer[100];
        n = sprintf(buffer, "%s[%i]", output_key.c_str(), int(jj));
        output_key = buffer;
      }

      PyObject *current_list = PyList_New(export_values[xx].size());
      for (size_t kk=0; kk<export_values[xx].size(); ++kk)
      {
        PyList_SetItem(current_list, kk, PyFloat_FromDouble(export_values[xx][kk]));
      }

      PyDict_SetItem(plot_dict, PyUnicode_FromString(output_key.c_str()), current_list);
      xx++;
    }
  }

  return plot_dict;
}


// Method definition
static PyMethodDef mod_methods[] = {{"parse_file", (PyCFunction)parse_file, METH_VARARGS|METH_KEYWORDS, "Method to parse a single silo file"},
                                {NULL, NULL, 0, NULL}};

// Module definition
static struct PyModuleDef mod_definition = {
    PyModuleDef_HEAD_INIT,
    "silo_parser",
    "Simple SILO parser",
    -1,
    mod_methods
};


PyMODINIT_FUNC PyInit_silo_parser(void)
{
  return PyModule_Create(&mod_definition);
}
