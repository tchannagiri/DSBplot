
SAM_MANDATORY_FIELDS = [
  'QNAME',
  'FLAG',
  'RNAME',
  'POS',
  'MAPQ',
  'CIGAR',
  'RNEXT',
  'PNEXT',
  'TLEN',
  'SEQ',
  'QUAL',
]

def parse_sam_fields(fields):
  """
  Get 2 dictionaries representing the information in the tab-separated SAM fields.

  Parameters
  ----------
  fields : a list of the tab-separated fields in one line of a SAM file.

  Returns
  -------
  mandatory : a dictionary of the mandatory SAM fields
    The dictionary keys are 'QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR',
    'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL'.
  optional :  a dictionary of the optional SAM fields
    Each dictionary key is the TAG ("AS", "XM", etc.,).
    Each value is itself a dictionary with keys "TYPE", "VALUE".
    Please refer to the SAM optional fields document (link below).
    
  Notes
  -----
  SAM file format : https://samtools.github.io/hts-specs/SAMv1.pdf
  SAM optional fields : https://github.com/samtools/hts-specs/blob/master/SAMtags.pdf
  """
  
  if len(fields) < len(SAM_MANDATORY_FIELDS):
    raise Exception('Not enough fields for SAM format')

  mandatory_fields = fields[:len(SAM_MANDATORY_FIELDS)]
  optional_fields = fields[len(SAM_MANDATORY_FIELDS):]
  mandatory = dict(zip(SAM_MANDATORY_FIELDS, mandatory_fields))
  optional = {}
  for field in optional_fields:
    field = field.split(':') # optional fields are of form TAG:TYPE:VALUE
    if len(field) != 3:
      raise Exception('Incorrect number of subfields: ' + str(field))
    optional[field[0]] = {'TYPE': field[1], 'VALUE': field[2]}
  return mandatory, optional