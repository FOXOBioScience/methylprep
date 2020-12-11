# Lib
import logging
import pandas as pd
# App
from ..models import (
    METHYLATED_PROBE_SUBSETS,
    UNMETHYLATED_PROBE_SUBSETS,
    METHYLATED_SNP_PROBES,
    UNMETHYLATED_SNP_PROBES,
    #METHYLATED_MOUSE_PROBES,
    #UNMETHYLATED_MOUSE_PROBES
)


__all__ = ['MethylationDataset']


LOGGER = logging.getLogger(__name__)


class MethylationDataset():
    """Wrapper for a collection of methylated or unmethylated probes and their mean intensity values,
    providing common functionality for the subset of probes.

    Arguments:
        raw_dataset {RawDataset} -- A sample's RawDataset for a single well on the processed array.
        manifest {Manifest} -- The Manifest for the correlated RawDataset's array type.
        probe_subsets {list(ProbeSubset)} -- Collection of ProbeSubsets that correspond to the probe type
        (methylated or unmethylated).
    """
    __bg_corrected = False
    __preprocessed = False # AKA NOOB CORRECTED

    def __init__(self, raw_dataset, manifest, probe_subsets):
        #LOGGER.info('Preprocessing methylation dataset: %s', raw_dataset.sample)

        self.probe_subsets = probe_subsets
        self.raw_dataset = raw_dataset # __init__ uses red_idat and green_idat IdatDatasets

        self.data_frames = {
            probe_subset: self._get_subset_means(manifest, probe_subset)
            for probe_subset in probe_subsets
        }

        self.data_frame = self.build_data_frame()

    @classmethod
    def methylated(cls, raw_dataset, manifest):
        """ convenience method that feeds in a pre-defined list of methylated CpG locii probes """
        return cls(raw_dataset, manifest, METHYLATED_PROBE_SUBSETS)

    @classmethod
    def unmethylated(cls, raw_dataset, manifest):
        """ convenience method that feeds in a pre-defined list of UNmethylated CpG locii probes """
        return cls(raw_dataset, manifest, UNMETHYLATED_PROBE_SUBSETS)

    @classmethod
    def snp_methylated(cls, raw_dataset, manifest):
        """ convenience method that feeds in a pre-defined list of methylated Snp locii probes """
        return cls(raw_dataset, manifest, METHYLATED_SNP_PROBES)

    @classmethod
    def snp_unmethylated(cls, raw_dataset, manifest):
        """ convenience method that feeds in a pre-defined list of UNmethylated Snp locii probes """
        return cls(raw_dataset, manifest, UNMETHYLATED_SNP_PROBES)

    #@classmethod
    #def mouse_methylated(cls, raw_dataset, manifest):
    #    """ convenience method that feeds in a pre-defined list of methylated MOUSE specific probes """
    #    return cls(raw_dataset, manifest, METHYLATED_MOUSE_PROBES)
    #
    #@classmethod
    #def mouse_unmethylated(cls, raw_dataset, manifest):
    #    """ convenience method that feeds in a pre-defined list of UNmethylated MOUSE specific probes """
    #    return cls(raw_dataset, manifest, UNMETHYLATED_MOUSE_PROBES)

    def build_data_frame(self):
        return pd.concat(self.data_frames.values())

    def _get_subset_means(self, manifest, probe_subset):
        channel_means = self.raw_dataset.get_channel_means(probe_subset.data_channel)
        channel_means = channel_means.assign(Channel=probe_subset.data_channel.value) # adds Red or Grn column

        probe_details = probe_subset.get_probe_details(manifest)
        #print(f"channel means: {channel_means.shape} -- ({probe_subset},{probe_subset.column_name}) -- probe_details {probe_details.shape}")
        # check here for probes that are missing data in manifest, and drop them if they are (better to be imperfect with warnings)
        #import pdb;pdb.set_trace()
        if probe_details[probe_subset.probe_address.header_name].isna().sum() > 0:
            print(f"probes missing from manifest: {probe_details[probe_subset.probe_address.header_name].isna().sum()}")
            print('These probes are probably incorrect in your manifest; processing cannot continue.')
            print( probe_details.loc[ probe_details[probe_subset.probe_address.header_name].isna() ].index )
            pre_shape = probe_details.shape
            probe_details = probe_details.drop( probe_details[ probe_details[probe_subset.probe_address.header_name].isna() ].index )
            print(f"{pre_shape[0] - probe_details.shape[0]} removed; {probe_details[probe_subset.probe_address.header_name].isna().sum()} nan remaining; but downstream steps will not work.")
            # this still won't fix it, because OOB also does some filtering.

        #print(f"probes in manifest missing here: ---TBD---")
        # -- no longer removes anything --- print(f"-- channel means: {channel_means.shape} -- REVISED probe_details {probe_details.shape}")

        # this merge matches AddressA_ID in probe_details to index (name: illumina_id) in channel_means.
        return probe_details.merge(
            channel_means,
            how='inner',
            left_on=probe_subset.probe_address.header_name,
            right_index=True,
            suffixes=(False, False),
        )

    def set_bg_corrected(self, green_corrected, red_corrected):
        # mouse array has duplicate readings for some 5000 probes; ignoring here.
        if green_corrected.duplicated(keep='first').sum() > 0:
            green_corrected.drop_duplicates(inplace=True)
        if red_corrected.duplicated(keep='first').sum() > 0:
            red_corrected.drop_duplicates(inplace=True)

        for probe_subset in self.data_frames:
            if probe_subset.is_red:
                corrected_values = red_corrected
            elif probe_subset.is_green:
                corrected_values = green_corrected
            else:
                raise ValueError('No data_channel for probe_subset')

            print(f"DEBUG NOOB shape mismatch: ({probe_subset},{probe_subset.column_name}) -- ({self.data_frames[probe_subset].shape} --> {corrected_values.loc[ self.data_frames[probe_subset][probe_subset.column_name] ].shape}) -vs- {corrected_values.shape}")
            #duplicate_values = green_corrected[ green_corrected['IlmnID'].isin( [k for k,v in dict(green_corrected['IlmnID'].value_counts()).items() if v>1] )]
            #PRE = self.data_frames[probe_subset].index
            #POST = corrected_values.loc[ self.data_frames[probe_subset][probe_subset.column_name] ].index
            #import pdb;pdb.set_trace()
            #print(f"-- ({probe_subset},{probe_subset.column_name}) PRE dupe {(set(PRE) - set(post))} POST dupe (set(POST) - set(PRE))")

            self._set_subset_bg_corrected(probe_subset, corrected_values)

        self.data_frame = self.build_data_frame()
        self.__bg_corrected = True

    def _set_subset_bg_corrected(self, probe_subset, corrected_values):
        # to begin, corrected_values has more rows, but only one copy of 'cg39970707_BC21'.
        # original [data_frame] has only one copy of 'cg39970707_BC21'
        # filtered_corrected has 5 copies of original.loc['cg39970707_BC21']
        original = self.data_frames[probe_subset]
        column = probe_subset.column_name
        filtered_corrected = corrected_values.loc[original[column]]

        #print(f"bg-correct {original.shape} [{probe_subset} {column}] --> {filtered_corrected.shape}")
        pre_filtered_corrected_shape = filtered_corrected.shape
        try:
            print(f"original cg39970707_BC21 {original.loc['cg39970707_BC21'].shape} --> corrected {filtered_corrected[filtered_corrected['IlmnID'] == 'cg39970707_BC21'].shape}")
        except KeyError:
            pass
        filtered_corrected = filtered_corrected.drop_duplicates()
        print(f"filtered_corrected {pre_filtered_corrected_shape} --deduped--> {filtered_corrected.shape}")

        updated = original.merge(
            filtered_corrected[['bg_corrected']],
            how='inner',
            left_on=column,
            right_index=True,
            suffixes=(False, False),
        )

        self.data_frames[probe_subset] = updated

    def set_noob(self, red_factor):
        for probe_subset, data_frame in self.data_frames.items():
            if probe_subset.is_red:
                data_frame = data_frame.assign(noob=data_frame['bg_corrected'] * red_factor)
            else:
                data_frame = data_frame.assign(noob=data_frame['bg_corrected'])

            self.data_frames[probe_subset] = data_frame

        self.data_frame = self.build_data_frame()
        self.__preprocessed = True
