""" BOADICEA web-service testing.  """
import os

from django.conf import settings
from django.contrib.auth.models import User
from django.urls import reverse
from django.test import TestCase
from rest_framework import status
from rest_framework.authtoken.models import Token
from rest_framework.test import APIClient
from django.utils.encoding import force_text
import json
from django.test.utils import override_settings
from datetime import date
from bws.pedigree import Female, CanRiskPedigree
from bws.cancer import Cancers, Cancer, CanRiskGeneticTests
from bws.calcs import Predictions


class MutFreqTests(TestCase):
    TEST_BASE_DIR = os.path.dirname(os.path.dirname(__file__))
    TEST_DATA_DIR = os.path.join(TEST_BASE_DIR, 'tests', 'data')

    @classmethod
    def setUpClass(cls):
        ''' Create a user, token and url. '''
        super(MutFreqTests, cls).setUpClass()

        cls.user = User.objects.create_user('testuser', email='testuser@test.com',
                                            password='testing')
        # add user details
        # UserDetails.objects.create(user=cls.user, job_title=UserDetails.CGEN,
        #                            country='UK')
        cls.user.save()
        cls.token = Token.objects.create(user=cls.user)
        cls.token.save()
        cls.url = reverse('bws')

        # create pedigree
        year = date.today().year
        gtests = CanRiskGeneticTests.default_factory()
        target = Female("FAM1", "F0", "001", "002", "003", target="1", age="42",
                        yob=str(year-42), cancers=Cancers(), gtests=gtests)
        cls.pedigree = CanRiskPedigree(people=[target])
        (_father, mother) = cls.pedigree.add_parents(target, gtests=gtests)
        mother.yob = str(year-67)
        mother.age = "67"
        mother.cancers = Cancers(bc1=Cancer("56"), bc2=Cancer(), oc=Cancer(), prc=Cancer(), pac=Cancer())

    def test_ashkn_mut_freq(self):
        '''
        Test POSTing CanRisk file with multiple families with and without Ashkenazi Jewish ancestry.
        Check mutation frequencies are correctly set in each case.
        '''
        client = APIClient(enforce_csrf_checks=True)
        pedigree_data = open(os.path.join(BwsTests.TEST_DATA_DIR, "canrisk_multi4xAJ.txt"), "r")

        data = {'mut_freq': 'UK', 'cancer_rates': 'UK',
                'pedigree_data': pedigree_data,
                'user_id': 'test_XXX'}
        client.credentials(HTTP_AUTHORIZATION='Token ' + MutFreqTests.token.key)
        response = client.post(MutFreqTests.url, data, format='multipart',
                               HTTP_ACCEPT="application/json")
        pedigree_data.close()
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        content = json.loads(force_text(response.content))
        self.assertTrue("mutation_frequency" in content)
        self.assertTrue("pedigree_result" in content)
        self.assertEqual(len(content['pedigree_result']), 4, "four results")

        family_ids = ["XXX2", "XX_AJ", "XX_NAJ", "XXX3"]
        for res in content['pedigree_result']:
            self.assertTrue(res['family_id'] in family_ids)

            if res['family_id'] == "XX_AJ" or res['family_id'] == "XX_NAJ":
                for r in res['cancer_risks']:
                    if r['age'] == 80:
                        c80 = r['breast cancer risk']['decimal']
                if res['family_id'] == "XX_AJ":
                    self.assertTrue("Ashkenazi" in res['mutation_frequency'])
                    c80_aj1 = c80
                elif res['family_id'] == "XX_NAJ":
                    self.assertTrue("UK" in res['mutation_frequency'])
                    c80_aj2 = c80

        self.assertGreater(c80_aj1, c80_aj2)
        self.assertTrue('mutation frequencies set to Ashkenazi Jewish population values for family (XX_AJ) as a ' +
                        'family member has Ashkenazi Jewish status.' in content['warnings'])

    @classmethod
    def get_risk(cls, crisks, age=80):
        ''' Return cancer risk to a specified age. '''
        for c in crisks:
            if c.get('age') == 80:
                return c
        return None

    @classmethod
    def get_mutation_prob(cls, mutation_probabilties, gene):
        ''' Return mutation probability for a specified gene. '''
        for mp in mutation_probabilties:
            this_gene = list(mp)[0]
            if this_gene == gene:
                return mp[this_gene]
        return None

    def test_custom_mutation_freq(self):
        ''' For different populations (UK, Canada), test CUSTOM mutation frequencies set to:
        1. same as default (UK) - check gives identical results
        2. different values to default - check gives different results
        '''
        ped = MutFreqTests.pedigree
        self.check_custom_mutations(settings.BC_MODEL, ped, 'breast cancer risk')
        self.check_custom_mutations(settings.OC_MODEL, ped, 'ovarian cancer risk')

    def check_custom_mutations(self, model, ped, ctext):
        '''
        @param model - cancer model settings (boadicea/ovarian)
        @param ped - Pedigree object
        @param ctext - cancer result text, e.g. breast cancer risk
        '''
        calcs = ['carrier_probs', 'remaining_lifetime']
        uk_mfreqs = model['MUTATION_FREQUENCIES']["UK"].copy()
        custom_mfreqs = model['MUTATION_FREQUENCIES']["UK"].copy()
        custom_mfreqs['BRCA1'] = 0.0012
        custom_mfreqs['BRCA2'] = 0.0015

        # cancer incidence rates
        cancer_rates = [model['CANCER_RATES'].get("UK"), model['CANCER_RATES'].get("Canada")]
        msens = model['GENETIC_TEST_SENSITIVITY']

        for cr in cancer_rates:
            calcs1 = Predictions(ped, cancer_rates=cr, calcs=calcs,     # 1. default frequencies (i.e. UK)
                                 mutation_sensitivity=msens, model_settings=model, mutation_frequency=uk_mfreqs)
            calcs2 = Predictions(ped, cancer_rates=cr, calcs=calcs,     # 2. custom frequencies same as default
                                 mutation_sensitivity=msens, model_settings=model, mutation_frequency=uk_mfreqs,
                                 population="Custom")
            calcs3 = Predictions(ped, cancer_rates=cr, calcs=calcs,     # 3. custom frequencies different to default
                                 mutation_sensitivity=msens, model_settings=model, population="Custom",
                                 mutation_frequency=custom_mfreqs)

            # lifetime risks to age 80
            bc1 = MutFreqTests.get_risk(calcs1.cancer_risks, age=80)[ctext]['decimal']
            bc2 = MutFreqTests.get_risk(calcs2.cancer_risks, age=80)[ctext]['decimal']
            bc3 = MutFreqTests.get_risk(calcs3.cancer_risks, age=80)[ctext]['decimal']

            self.assertEqual(bc1, bc2, "Lifetime risk to 80 same")
            self.assertNotEqual(bc1, bc3, "Lifetime risk to 80 different")

            # gene mutation probabilities
            mprob = MutFreqTests.get_mutation_prob
            for gene in model['GENES']:
                self.assertEqual(mprob(calcs1.mutation_probabilties, gene)['decimal'],
                                 mprob(calcs2.mutation_probabilties, gene)['decimal'])

            self.assertNotEqual(mprob(calcs1.mutation_probabilties, 'BRCA1')['decimal'],
                                mprob(calcs3.mutation_probabilties, 'BRCA1')['decimal'])
            self.assertNotEqual(mprob(calcs1.mutation_probabilties, 'BRCA2')['decimal'],
                                mprob(calcs3.mutation_probabilties, 'BRCA2')['decimal'])

    def test_custom_mutation_freq_iceland(self):
        ''' For Iceland populations, test CUSTOM mutation frequencies set to:
        1. same as default (Iceland) - check gives identical results
        2. different values to default - check gives different results
        '''
        calcs = ['remaining_lifetime']
        ped = MutFreqTests.pedigree
        models = [settings.BC_MODEL, settings.OC_MODEL]     # boadicea and ovarian models
        for model in models:
            ice_mfreqs = model['MUTATION_FREQUENCIES']["Iceland"].copy()
            custom_mfreqs = model['MUTATION_FREQUENCIES']["Iceland"].copy()
            custom_mfreqs['BRCA1'] = 0.0012

            ice_cr = model['CANCER_RATES'].get("Iceland")

            msens = model['GENETIC_TEST_SENSITIVITY']

            calcs1 = Predictions(ped, cancer_rates=ice_cr, calcs=calcs,     # 1. default frequencies (i.e. Iceland)
                                 mutation_frequency=ice_mfreqs,
                                 mutation_sensitivity=msens, model_settings=model)
            calcs2 = Predictions(ped, cancer_rates=ice_cr, calcs=calcs,     # 2. custom frequencies same as default
                                 mutation_frequency=ice_mfreqs, population="Custom",
                                 mutation_sensitivity=msens, model_settings=model)
            calcs3 = Predictions(ped, cancer_rates=ice_cr, calcs=calcs,     # 3. custom frequencies different to default
                                 mutation_frequency=custom_mfreqs, population="Custom",
                                 mutation_sensitivity=msens, model_settings=model)

            ctext = 'breast cancer risk' if model['NAME'] == 'BC' else 'ovarian cancer risk'
            c1 = MutFreqTests.get_risk(calcs1.cancer_risks, age=80)[ctext]['decimal']
            c2 = MutFreqTests.get_risk(calcs2.cancer_risks, age=80)[ctext]['decimal']
            c3 = MutFreqTests.get_risk(calcs3.cancer_risks, age=80)[ctext]['decimal']
            self.assertEqual(c1, c2, "Lifetime risk to 80 same")
            self.assertNotEqual(c1, c3, "Lifetime risk to 80 different")


class BwsTests(TestCase):
    TEST_BASE_DIR = os.path.dirname(os.path.dirname(__file__))
    TEST_DATA_DIR = os.path.join(TEST_BASE_DIR, 'tests', 'data')

    @classmethod
    def setUpClass(cls):
        ''' Create a user, token and url. '''
        super(BwsTests, cls).setUpClass()

        cls.user = User.objects.create_user('testuser', email='testuser@test.com',
                                            password='testing')
        # add user details
        # UserDetails.objects.create(user=cls.user, job_title=UserDetails.CGEN,
        #                            country='UK')
        cls.user.save()
        cls.token = Token.objects.create(user=cls.user)
        cls.token.save()
        cls.url = reverse('bws')

    def setUp(self):
        ''' Set up test client and pedigree data. '''
        self.client = APIClient(enforce_csrf_checks=True)
        self.pedigree_data = open(os.path.join(BwsTests.TEST_DATA_DIR, "pedigree_data.txt"), "r")

    def tearDown(self):
        TestCase.tearDown(self)
        self.pedigree_data.close()

    def test_token_auth_bws(self):
        ''' Test POSTing to the BWS using token authentication. '''
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK',
                'pedigree_data': self.pedigree_data,
                'user_id': 'test_XXX'}
        self.client.credentials(HTTP_AUTHORIZATION='Token ' + BwsTests.token.key)
        response = self.client.post(BwsTests.url, data, format='multipart',
                                    HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        content = json.loads(force_text(response.content))
        self.assertTrue("mutation_frequency" in content)
        self.assertTrue("pedigree_result" in content)
        self.assertTrue("family_id" in content["pedigree_result"][0])

    def test_multi_pedigree_bws(self):
        ''' Test POSTing multiple pedigrees to the BWS. '''
        multi_pedigree_data = open(os.path.join(BwsTests.TEST_DATA_DIR, "multi_pedigree_data.txt"), "r")
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK',
                'pedigree_data': multi_pedigree_data,
                'user_id': 'test_XXX'}
        self.client.credentials(HTTP_AUTHORIZATION='Token ' + BwsTests.token.key)
        response = self.client.post(BwsTests.url, data, format='multipart',
                                    HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        content = json.loads(force_text(response.content))
        self.assertEqual(len(content['pedigree_result']), 2, "two results")
        family_ids = ["XXX0", "XXX1"]
        for res in content['pedigree_result']:
            self.assertTrue(res['family_id'] in family_ids)
        multi_pedigree_data.close()

    def test_canrisk_format_bws(self):
        ''' Test POSTing canrisk format pedigree to the BWS. '''
        canrisk_data = open(os.path.join(BwsTests.TEST_DATA_DIR, "canrisk_data_v1.txt"), "r")
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK',
                'pedigree_data': canrisk_data,
                'user_id': 'test_XXX'}
        self.client.credentials(HTTP_AUTHORIZATION='Token ' + BwsTests.token.key)
        response = self.client.post(BwsTests.url, data, format='multipart',
                                    HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        content = json.loads(force_text(response.content))
        self.assertTrue("pedigree_result" in content)
        genes = settings.BC_MODEL['GENES']
        for g in genes:
            self.assertTrue(g in content['mutation_frequency']['UK'])
        canrisk_data.close()

    def test_token_auth_err(self):
        ''' Test POSTing to the BWS using token authentication. '''
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK',
                'pedigree_data': self.pedigree_data}
        self.client.credentials(HTTP_AUTHORIZATION='Token xxxxxxxxxx')
        response = self.client.post(BwsTests.url, data, format='multipart',
                                    HTTP_ACCEPT="application/json")
        self.assertNotEqual(response.status_code, status.HTTP_200_OK)
        content = json.loads(force_text(response.content))
        self.assertEqual(content['detail'], 'Invalid token.')

    def test_force_auth_bws(self):
        ''' Test POSTing to the BWS bypassing authentication. '''
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK',
                'pedigree_data': self.pedigree_data,
                'user_id': 'test_XXX'}
        self.client.force_authenticate(user=BwsTests.user)
        response = self.client.post(BwsTests.url, data, format='multipart')
        self.assertEqual(response.status_code, status.HTTP_200_OK)

    def test_custom_mutation_frequency(self):
        ''' Test POSTing custom mutation frequencies. '''
        genes = settings.BC_MODEL['GENES']
        data = {g.lower() + '_mut_frequency': 0.00085 for g in genes}
        data.update({'mut_freq': 'Custom', 'cancer_rates': 'UK',
                     'pedigree_data': self.pedigree_data, 'user_id': 'test_XXX'})

        self.client.force_authenticate(user=BwsTests.user)
        response = self.client.post(BwsTests.url, data, format='multipart',
                                    HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        content = json.loads(force_text(response.content))

        for g, mf in content['mutation_frequency']['Custom'].items():
            self.assertTrue(g in genes)
            self.assertEqual(mf, data[g.lower() + '_mut_frequency'])

    def test_custom_mutation_frequency_errs(self):
        ''' Test POSTing custom mutation frequencies with errors. '''
        genes = settings.BC_MODEL['GENES']
        data = {g.lower() + '_mut_frequency':
                (settings.MAX_MUTATION_FREQ + 0.1) if idx % 2 == 0 else (settings.MAX_MUTATION_FREQ - 0.1)
                for idx, g in enumerate(genes)}
        data.update({'mut_freq': 'Custom', 'cancer_rates': 'UK',
                     'pedigree_data': self.pedigree_data, 'user_id': 'test_XXX'})

        self.client.force_authenticate(user=BwsTests.user)
        response = self.client.post(BwsTests.url, data, format='multipart',
                                    HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
        content = json.loads(force_text(response.content))
        self.assertEqual(len(content.keys()), len(genes))
        for k in content.keys():
            self.assertTrue(k.split("_")[0].upper() in genes)

    def test_missing_fields(self):
        ''' Test POSTing with missing required fields. '''
        data = {'mut_freq': 'UK'}
        self.client.credentials(HTTP_AUTHORIZATION='Token ' + BwsTests.token.key)
        response = self.client.post(BwsTests.url, data, format='multipart',
                                    HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
        content = json.loads(force_text(response.content))
        self.assertEqual(content['user_id'][0], 'This field is required.')
        self.assertEqual(content['cancer_rates'][0], 'This field is required.')
        self.assertEqual(content['pedigree_data'][0], 'This field is required.')

    def test_bws_errors(self):
        ''' Test an error is reported by the web-service for an invalid year of birth. '''
        # force an error changing to an invalid year of birth
        ped = open(os.path.join(BwsTests.TEST_DATA_DIR, "pedigree_data.txt"), "r")
        pd = ped.read().replace('1963', '1600')
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK',
                'pedigree_data': pd, 'user_id': 'test_XXX'}
        self.client.force_authenticate(user=BwsTests.user)
        response = self.client.post(BwsTests.url, data, format='multipart',
                                    HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
        content = json.loads(force_text(response.content))
        ped.close()
        self.assertTrue('Person Error' in content)
        self.assertTrue('year of birth' in content['Person Error'])

    @override_settings(FORTRAN_TIMEOUT=0.05)
    def test_bws_timeout(self):
        ''' Test a timeout error is reported by the web-service. '''
        data = {'mut_freq': 'UK', 'cancer_rates': 'UK',
                'pedigree_data': self.pedigree_data, 'user_id': 'test_XXX'}
        self.client.force_authenticate(user=BwsTests.user)
        response = self.client.post(BwsTests.url, data, format='multipart',
                                    HTTP_ACCEPT="application/json")
        self.assertEqual(response.status_code, status.HTTP_408_REQUEST_TIMEOUT)
        content = json.loads(force_text(response.content))
        self.assertTrue('detail' in content)
        self.assertTrue('Request has timed out.' in content['detail'])
