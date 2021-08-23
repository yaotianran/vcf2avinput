#!/usr/bin/env python3
# Used in somatic mutation calling pipeline.
# Remove the records without PASS and tumor support reads less than 3 float
# 0.5
import os
import sys
import argparse
import os.path as path
import scipy.stats as stats
import vcfpy

from math import *

parser_ar = argparse.ArgumentParser(prog='vcf2avinput.py',
                                    usage='vcf2avinput.py [OPTION] <vcf_file>',
                                    description='To replace the convert2annovar.pl from ANNOVAR. Multiple samples supported.\nUse bcftools to normalize your vcf file first.',
                                    epilog='When use table_annovar.pl, add --otherinfo option to retain the additional information\n',
                                    formatter_class=argparse.RawTextHelpFormatter)

parser_ar.add_argument('vcf_file', help='FILE. A vcf file.', metavar='vcf_file')

parser_ar.add_argument('-t', metavar='auto', default='auto', choices = ['auto', 'hap', 'mutect', 'ss', 'si', 'deepvariant', 'varscan'], help='STR. The type of VCF file. Currently can be assigned as one of [hap, mutect, ss, si, deepvariant, varscan].', dest='TYPE')
parser_ar.add_argument('-o', default='', help='Output avinput file name.', dest='OUTPUT')
parser_ar.add_argument('-s', metavar='"sample_name_1,..."', default='', help='Sample name. Multiple samples can be seperated by comma. If not assigned, script will extract every samples in the vcf file.', dest='SAMPLE_NAME')

parser_ar.add_argument('-q', metavar='[20.0]', default=20.0, type=float, help='FLOAT. QUAL cutoff.', dest='QUAL_CUTOFF')
parser_ar.add_argument('-d', metavar='[0]', default=0, type=int, help='INT. The read depth cutoff, equals to DP in FORMAT column or AN in INFO column.', dest='DEPTH')
parser_ar.add_argument('-c', metavar='[0]', default=0, type=int, help='INT. The alt allele supported read depth, equals to AD in FORMAT column or AC in INFO column.', dest='ALLELE_DEPTH')
parser_ar.add_argument('-f', metavar='[0.0]', default=0.0, type=float, help='FLOAT. The minimum allele frequence. Calculate from AD/DP or AF in INFO column.', dest='ALLELE_FREQUENCE')
parser_ar.add_argument('-pr', metavar='[0]', default=0, type=int, help='INT. The pair read count supporting an alternative SV allele.', dest='PAIR_READ')
parser_ar.add_argument('-sr', metavar='[0]', default=0, type=int, help='INT. The split read count supporting an alternative SV allele.', dest='SPLIT_READ')

parser_ar.add_argument('--all', action='store_true', default=False, help='By defalut the script will only output PASS variants, turning on this switch will let you output all variants', dest='OUTPUT_ALL')
parser_ar.add_argument('--fail_on_missing', action='store_true', default=False, help='By defalut the script treat variants with missing filters as PASS, the switch will treat them as FAIL.', dest='FAIL_ON_MISSING')
parser_ar.add_argument('--remove_ao', action='store_true', default=False, help='For lumpy type only, remove any call with AO (Alternate allele observations)', dest='REMOVE_AO')

args = parser_ar.parse_args()

VCFFILE = args.vcf_file
TYPE = args.TYPE
SAMPLE_NAME_LST = args.SAMPLE_NAME.split(',')
OUTPUT_FILE = args.OUTPUT

QUAL_CUTOFF = args.QUAL_CUTOFF
DEPTH = args.DEPTH
ALLELE_DEPTH = args.ALLELE_DEPTH
ALLELE_FREQUENCE = args.ALLELE_FREQUENCE
PAIR_READ = args.PAIR_READ
SPLIT_READ = args.SPLIT_READ

OUTPUT_ALL = args.OUTPUT_ALL
FAIL_ON_MISSING = args.FAIL_ON_MISSING
REMOVE_AO = args.REMOVE_AO


class vcf_record:

   record = None  # vcfpy.record
   #__chrom_str = ''
   #__start_int = 0
   #__end_int = 0
   #__ref_str = ''
   #__alt_str = ''

   def __init__(self, record):
      self.record = record
      if len(self.record.ALT) > 1:
         message = 'The input VCF file can not include multi-allelics records.\n{}:{} {} -> {}'.format(self.record.CHROM,
                                                                                                            self.record.POS,
                                                                                                            self.record.REF,
                                                                                                            self.record.ALT)
         raise ValueError(message)

      if self.record.ALT[0].type not in ['SNV', 'MNV', 'INS', 'DEL', 'INDEL']:  # ['BND', 'SYMBOLIC']:
         message = 'This record contain unknown types.\n{}:{} {} -> {}'.format(self.record.CHROM,
                                                                               self.record.POS,
                                                                               self.record.REF,
                                                                               self.record.ALT)
         raise TypeError(message)

      return None


   @property
   def type(self) -> str:
      return self.record.ALT[0].type

   @property
   def chrom(self) -> str:
      return self.record.CHROM

   @property
   def start(self) -> int:
      start = 0
      if self.type == 'SNV' or self.type == 'MNV' or self.type == 'INS' or self.type == 'INDEL':
         start = self.record.POS
      elif self.type == 'DEL':
         start = self.record.POS + 1

      return start

   @property
   def end(self):
      end = 0
      if self.type == 'SNV' or self.type == 'MNV' or self.type == 'INS' or self.type == 'DEL' or self.type == 'INDEL':
         end = self.record.POS + len(self.record.REF) - 1

      return end

   @property
   def ref(self) -> str:
      ref = ''
      if self.type == 'SNV' or self.type == 'MNV' or self.type == 'INDEL':
         ref = self.record.REF
      elif self.type == 'INS':
         ref = '-'
      elif self.type == 'DEL':
         ref = self.record.REF[1:]
      else:
         pass

      return ref

   @property
   def alt(self) -> str:
      alt = 'ini'
      if self.type == 'SNV' or self.type == 'MNV' or self.type == 'INDEL':
         alt = self.record.ALT[0].value
      elif self.type == 'INS':
         alt = self.record.ALT[0].value[1:]
      elif self.type == 'DEL':
         alt = '-'
      else:
         pass

      return alt

   @property
   def qual(self) -> int:
      return self.record.QUAL

   @property
   def filter(self) -> list:
      return self.record.FILTER


class vcf_record_HaplotypeCaller(vcf_record):

   # {'sample_name', [AD:int, DP:int, GT:str:{het, hom}, GQ:int]}
   @property
   def GT(self) -> dict:

      sample_name = ''
      AD = 0
      DP = 0
      GT = ''
      GQ = 0

      result_dict = {}
      for call in self.record.calls:
         if call.data['GT'] in ['0/1', '1/1']:
            sample_name = call.sample

            try:
               AD = call.data['AD'][1]
            except:
               AD = None

            try:
               DP = call.data['DP']
            except:
               DP = None

            GT = 'hom' if call.data['GT'] == '1/1' else 'het'

            try:
               GQ = call.data['GQ']
            except:
               GQ = None

            result_dict[sample_name] = [AD, DP, GT, GQ]

      return result_dict

class vcf_record_mutect2(vcf_record):

   # {'sample_name', [AD:int, DP:int, GT:str:{het, hom}, weighted_average_AF:float, FR_bias:float]}
   # FR_bias (0-1) should not larger than 0.8
   @property
   def GT(self) -> dict:
      sample_name = ''
      AD = 0
      DP = 0
      GT = ''
      # GQ = 0
      weighted_average_AF = 0.0
      FR_bias = 0.0

      result_dict = {}
      for call in self.record.calls:
         if call.data['GT'] in ['0/1', '1/1']:
            sample_name = call.sample

            try:
               AD = call.data['AD'][1]
            except:
               AD = None

            try:
               DP = call.data['DP']
            except:
               DP = None

            GT = 'hom' if call.data['GT'] == '1/1' else 'het'

            #try:
               #GQ = call.data['GQ']
            #except:
               #GQ = None

            try:
               weighted_average_AF = call.data['AF'][0]
            except:
               weighted_average_AF = None

            try:
               # FR_bias = float(format(sqrt((call.data['F1R2'][1] - AD / 2) ** 2 + (call.data['F2R1'][1] - AD / 2) ** 2) / AD, '.3g'))
               FR_bias =  abs(call.data['F1R2'][1] - call.data['F2R1'][1]) / (call.data['F1R2'][1] + call.data['F2R1'][1])
               FR_bias = float(format(FR_bias, '.3g'))
            except:
               FR_bias = None

            result_dict[sample_name] = [AD, DP, GT, weighted_average_AF, FR_bias]

      return result_dict


class vcf_record_strelka_snv(vcf_record):

   # {'sample_name', [AD:int, DP:int, GT:str:{het, hom}, QSS:int, PVAL:float]}
   # SNVSB = Somatic SNV site strand bias (should not larger than 10)
   # QSS = Quality score for any somatic snv
   # PVAL = p-value for Fisher exact test
   @property
   def GT(self) -> dict:
      sample_name = ''
      AD = 0
      DP = 0
      GT = ''

      QSS = 0
      PVAL = 0.0

      result_dict = {}

      try:
         GT = self.record.INFO['SGT'].split('>')[1]
         GT = 'hom' if GT in ['AA', 'TT', 'CC', 'GG'] else 'het'
      except:
         GT = None

      #try:
         #SNVSB = self.record.INFO['SNVSB']
      #except:
         #SNVSB = None

      try:
         QSS = self.record.INFO['QSS']
      except:
         QSS = None

      a = 0
      b = 0
      c = 0
      d = 0
      for call in self.record.calls:
         sample_name = call.sample

         try:
            if self.record.ALT[0].value == 'A':
               AD = call.data['AU'][0]
            elif self.record.ALT[0].value == 'T':
               AD = call.data['TU'][0]
            elif self.record.ALT[0].value == 'C':
               AD = call.data['CU'][0]
            elif self.record.ALT[0].value == 'G':
               AD = call.data['GU'][0]
            else:
               raise ValueError
         except:
            AD = None

         try:
            DP = call.data['DP'] - call.data['FDP']
         except:
            DP = None

         result_dict[sample_name] = [AD, DP, GT, QSS]

         try:
            if sample_name == 'NORMAL':
               a = AD
               b = DP
            elif sample_name == 'TUMOR':
               c = AD
               d = DP
            else:
               raise ValueError
         except:
            pass

      try:
         r = stats.fisher_exact([[a, b], [c, d]], alternative = 'less')
         PVAL = float(format(-log10(r[1]), '.3g'))
      except:
         PVAL = None

      for key in result_dict.keys():
         result_dict[key].append(PVAL)

      return result_dict

class vcf_record_strelka_indel(vcf_record):

   # {'sample_name', [AD:int, DP:int, GT:str:{het, hom}, QSI:int, PVAL:float]}
   @property
   def GT(self) -> dict:
      sample_name = ''
      AD = 0
      DP = 0
      GT = ''

      QSI = 0
      PVAL = 0.0

      result_dict = {}

      try:
         GT = self.record.INFO['SGT'].split('>')[1]
      except:
         GT = None

      try:
         QSI = self.record.INFO['QSI']
      except:
         QSI = None

      a = 0
      b = 0
      c = 0
      d = 0
      for call in self.record.calls:
         sample_name = call.sample

         try:
            AD = call.data['TIR'][0]
         except:
            AD = None

         try:
            DP = call.data['TIR'][0] + call.data['TAR'][0]
         except:
            DP = None

         result_dict[sample_name] = [AD, DP, GT, QSI]

         try:
            if sample_name == 'NORMAL':
               a = AD
               b = DP
            elif sample_name == 'TUMOR':
               c = AD
               d = DP
            else:
               raise ValueError
         except:
            pass

      try:
         r = stats.fisher_exact([[a, b], [c, d]], alternative = 'less')
         PVAL = float(format(-log10(r[1]), '.3g'))
      except:
         PVAL = None

      for key in result_dict.keys():
         result_dict[key].append(PVAL)

      return result_dict


class vcf_record_deepvariant(vcf_record):

   # {'sample_name', [AD:int, DP:int, GT:str:{het, hom}, GQ:int]}
   @property
   def GT(self) -> dict:

      sample_name = ''
      AD = 0
      DP = 0
      GT = ''
      GQ = 0

      result_dict = {}
      for call in self.record.calls:
         if call.data['GT'] in ['0/1', '1/1']:
            sample_name = call.sample

            try:
               AD = call.data['AD'][1]
            except:
               AD = None

            try:
               DP = call.data['DP']
            except:
               DP = None

            GT = 'hom' if call.data['GT'] == '1/1' else 'het'

            try:
               GQ = call.data['GQ']
            except:
               GQ = None

            result_dict[sample_name] = [AD, DP, GT, GQ]

      return result_dict

class vcf_record_varscan(vcf_record):

   # {'sample_name', [AD:int, DP:int, GT:str:{het, hom}, GQ:int, FR_bias:float, PVAL:str]}
   @property
   def GT(self) -> dict:

      sample_name = ''
      AD = 0
      DP = 0
      GT = ''
      GQ = 0
      FR_bias = 0.0
      p_value = 0

      result_dict = {}
      for call in self.record.calls:
         if call.data['GT'] in ['0/1', '1/1']:
            sample_name = call.sample

            try:
               AD = call.data['AD'][0]
            except:
               AD = None

            try:
               DP = call.data['DP']
            except:
               DP = None

            GT = 'hom' if call.data['GT'] == '1/1' else 'het'

            try:
               GQ = call.data['GQ']
            except:
               GQ = None

            try:
               # FR_bias = float(format(sqrt((call.data['F1R2'][1] - AD / 2) ** 2 + (call.data['F2R1'][1] - AD / 2) ** 2) / AD, '.3g'))
               FR_bias = abs(call.data['ADF'] - call.data['ADR']) / (call.data['ADF'] + call.data['ADR'])
               FR_bias = float(format(FR_bias, '.3g'))
            except:
               FR_bias = None

            try:
               p_value = float(format(-log10(float(call.data['PVAL'])), '.3g'))
               p_value = min([p_value, 255])
            except:
               p_value = None

            result_dict[sample_name] = [AD, DP, GT, GQ, FR_bias, p_value]

      return result_dict



def detect_type(vcf_file: str) -> str:
   # hap, mutect, ss, si, deepvariant, varscan
   type_str = 'unknown'

   vcf = vcfpy.Reader.from_path(path.realpath(path.expanduser(VCFFILE)))
   for line in vcf.header.lines:  # line = vcf.HeadLine
      try:
         if line.key == 'source':
            if line.value == 'HaplotypeCaller':
               return 'hap'
            if line.value == 'Mutect2':
               return 'mutect'
            if line.value == 'strelka':
               if 'QSS' in vcf.header.info_ids():
                  return 'ss'
               if 'QSI' in vcf.header.info_ids():
                  return 'si'
            if line.value == 'VarScan2':
               return 'varscan'
         elif line.key == 'DeepVariant_version':
            return 'deepvariant'
         elif line.key == 'MutectVersion':
            return 'mutect'
         else:
            continue


      except:
         continue
   return type_str

def main(argvList=sys.argv, argv_int=len(sys.argv)):

   # get input
   try:
      vcf = vcfpy.Reader.from_path(path.realpath(path.expanduser(VCFFILE)))
   except Exception as ex:
      print(ex)
      sys.exit(1)

   # get output
   if OUTPUT_FILE == '':
      output_file = path.splitext(path.realpath(path.expanduser(VCFFILE)))[0] + '.avinput'
   else:
      output_file = path.realpath(path.expanduser(OUTPUT_FILE))

   # get sample
   sample_lst = []
   if SAMPLE_NAME_LST == ['']:
      sample_lst = vcf.header.samples.names
   else:
      sample_lst = SAMPLE_NAME_LST
   print('{} has {} samples: {}'.format(VCFFILE, len(vcf.header.samples.names), vcf.header.samples.names))
   print('Extracting samples: {}'.format(sample_lst))

   # ================================================================
   file_type = ''
   if TYPE == 'auto':  # hap, mutect, ss, si, deepvariant, varscan
      file_type = detect_type(VCFFILE)
      if file_type == 'unknown':
         message = 'Can not determine the input vcf type.\nPlease use -t parameter to specify a type [hap, mutect, ss, si, deepvariant, varscan].'
         raise TypeError(message)
      else:
         message = 'Detect the vcf file as {} type.'.format(file_type)
         print(message)
   else:
      file_type = TYPE
      message = 'Treat the vcf file as {} type.'.format(file_type)
      print(message)

   with open(output_file, 'wt', buffering=1) as output_f:
      i = 0
      for record in vcf:
         i += 1
         if i % 1000 == 0:
            print(i, end = '\r')

         # hap, mutect, ss, si, deepvariant, varscan
         if file_type == 'hap':
            r = vcf_record_HaplotypeCaller(record)
         elif file_type == 'mutect':
            r = vcf_record_mutect2(record)
         elif file_type == 'ss':
            r = vcf_record_strelka_snv(record)
         elif file_type == 'si':
            r = vcf_record_strelka_indel(record)
         elif file_type == 'deepvariant':
            r = vcf_record_deepvariant(record)
         elif file_type == 'varscan':
            r = vcf_record_varscan(record)
         else:
            message = ''
            raise TypeError(message)

         # ==============================
         if r.filter == [] and FAIL_ON_MISSING:
            continue

         if r.filter != ['PASS'] and not OUTPUT_ALL:
            continue

         if isinstance(r.qual, (int, float)) and r.qual < QUAL_CUTOFF:
            continue

         if r.alt == '*':
            message = 'The allele for this record is missed due to a upstream deletion.\n{}:{} {} -> {}'.format(r.record.CHROM,
                                                                                                                  r.record.POS,
                                                                                                                  r.record.REF,
                                                                                                                  r.record.ALT)
            print(message)
            continue

         line_str = ''
         for key, value in r.GT.items():  # key is a sample name, value is a list
            if key not in sample_lst:
               continue

            info_str = '\t'.join(map(str, value))
            line_str = '{chrom}\t{start}\t{end}\t{ref}\t{alt}\t{sample}\t{info}'.format(chrom = r.chrom,
                                                                                        start = r.start,
                                                                                        end = r.end,
                                                                                        ref = r.ref,
                                                                                        alt = r.alt,
                                                                                        sample = key,
                                                                                        info = info_str)
            output_f.writelines(line_str + '\n')
      print('{} records proceeded and save to {}.'.format(i, output_file))
      print('Done')

   return


if __name__ == '__main__':
   main()


