import json
import time
from datetime import datetime
from typing import Dict, List, Optional, Set, Tuple
from itertools import combinations
from dataclasses import dataclass, asdict
import pytest
from agent_management.agents.pubchem_agent import PubChemAgent
from agent_management.llm_service import LLMService, LLMModelConfig, ProviderType
import os

@dataclass
class MoleculeResult:
    name: str
    cid: int
    formula: str

@dataclass
class SearchMethodResult:
    success: bool
    interpreted_molecule: str
    num_molecules_found: int
    molecules: List[MoleculeResult]
    error: Optional[str]
    processing_time: float

@dataclass
class TopicResults:
    topic: str
    direct_rest_exact: SearchMethodResult
    direct_rest_fuzzy: SearchMethodResult
    pubchempy_name: SearchMethodResult
    pubchempy_formula: SearchMethodResult
    sourcetable: SearchMethodResult

@dataclass
class TestResults:
    total_topics: int
    timestamp: str
    topics: List[TopicResults]
    
    def calculate_stats(self) -> Dict:
        stats = {}
        methods = ['direct_rest_exact', 'direct_rest_fuzzy', 'pubchempy_name', 'pubchempy_formula', 'sourcetable']
        
        for method in methods:
            successes = sum(1 for topic in self.topics if getattr(topic, method).success)
            stats[f"{method}_success_rate"] = (successes / self.total_topics) * 100
            stats[f"{method}_avg_time"] = sum(getattr(topic, method).processing_time for topic in self.topics) / self.total_topics
            
            # Calculate unique successes
            unique_successes = sum(1 for topic in self.topics 
                if getattr(topic, method).success and 
                sum(1 for other_method in methods if other_method != method and getattr(topic, other_method).success) == 0)
            stats[f"{method}_unique_successes"] = unique_successes
        
        # Calculate best method combinations
        for i in range(2, len(methods) + 1):
            for combo in combinations(methods, i):
                combo_successes = sum(1 for topic in self.topics 
                    if any(getattr(topic, method).success for method in combo))
                stats[f"{'_'.join(combo)}_combined_success_rate"] = (combo_successes / self.total_topics) * 100
        
        return stats

# Test topics with their expected molecule names
TEST_TOPICS = [
    ("norbornane", "Teach me about norbornane and its bridging structure"),
    ("cryptand-222", "Teach me about cryptand[2.2.2] and 3D host-guest complexation"),
    ("Mo2(O2CCH3)4", "Teach me about metal-metal quadruple bonds in dimolybdenum complexes"),
    ("closo-B12H12-2", "Teach me about carboranes and their polyhedral cage structures"),
    ("[2]rotaxane", "Teach me about rotaxanes and their mechanically interlocked architecture"),
    ("PAMAM", "Teach me about dendrimers and their hyperbranched growth patterns"),
    ("K3C60", "Teach me about alkali-doped fullerenes and superconductivity"),
    ("adenine", "Teach me about the DNA double helix and base-pair stacking"),
    ("copper phthalocyanine", "Teach me about metallophthalocyanines and their planar macrocycles"),
    ("phosphotungstate", "Teach me about polyoxometalates and their metal-oxygen clusters")
]

def process_topic_with_method(agent: PubChemAgent, molecule_name: str, topic: str, method: str, **kwargs) -> SearchMethodResult:
    start_time = time.time()
    try:
        # Use the agent's normalize_query method to get variations
        query_variations = agent._normalize_query(molecule_name)
        
        # Try each variation until we get a result
        for query in query_variations:
            try:
                if method == 'direct_rest_exact':
                    molecules = agent._search_pubchem_direct(query)
                    if molecules:
                        break
                elif method == 'direct_rest_fuzzy':
                    molecules = agent._search_pubchem_rest(query)
                    if molecules:
                        break
                elif method == 'pubchempy_name':
                    molecules = agent._search_pubchem_direct(query)
                    if molecules:
                        break
                elif method == 'pubchempy_formula':
                    molecules = agent._search_pubchem_rest(query)
                    if molecules:
                        break
                elif method == 'sourcetable':
                    molecules = agent._search_pubchem_direct(query)
                    if molecules:
                        break
                else:
                    raise ValueError(f"Unknown method: {method}")
            except Exception as e:
                molecules = []
                continue
        
        molecule_results = []
        for mol in molecules:
            # Handle both dictionary and Compound objects
            if isinstance(mol, dict):
                molecule_results.append(MoleculeResult(
                    name=mol.get('name', ''),
                    cid=mol.get('cid', 0),
                    formula=mol.get('formula', '')
                ))
            else:
                molecule_results.append(MoleculeResult(
                    name=mol.name if hasattr(mol, 'name') else '',
                    cid=mol.cid if hasattr(mol, 'cid') else 0,
                    formula=mol.molecular_formula if hasattr(mol, 'molecular_formula') else ''
                ))
        
        return SearchMethodResult(
            success=len(molecule_results) > 0,
            interpreted_molecule=molecule_name,
            num_molecules_found=len(molecule_results),
            molecules=molecule_results,
            error=None,
            processing_time=time.time() - start_time
        )
    except Exception as e:
        return SearchMethodResult(
            success=False,
            interpreted_molecule=molecule_name,
            num_molecules_found=0,
            molecules=[],
            error=str(e),
            processing_time=time.time() - start_time
        )

def test_search_method_comparison():
    # Initialize LLM service with default configuration
    config = LLMModelConfig(
        provider=ProviderType.OPENAI,  # Default to OpenAI
        model_name="gpt-3.5-turbo",  # Default model
        api_key=None  # Will use environment variable
    )
    llm_service = LLMService(config)
    agent = PubChemAgent(llm_service)
    
    results = TestResults(
        total_topics=len(TEST_TOPICS),
        timestamp=datetime.now().isoformat(),
        topics=[]
    )
    
    for molecule_name, topic in TEST_TOPICS:
        topic_results = TopicResults(
            topic=topic,
            direct_rest_exact=process_topic_with_method(agent, molecule_name, topic, 'direct_rest_exact'),
            direct_rest_fuzzy=process_topic_with_method(agent, molecule_name, topic, 'direct_rest_fuzzy'),
            pubchempy_name=process_topic_with_method(agent, molecule_name, topic, 'pubchempy_name'),
            pubchempy_formula=process_topic_with_method(agent, molecule_name, topic, 'pubchempy_formula'),
            sourcetable=process_topic_with_method(agent, molecule_name, topic, 'sourcetable')
        )
        results.topics.append(topic_results)
    
    # Calculate statistics
    stats = results.calculate_stats()
    
    # Save detailed results
    with open('search_methods_comparison.json', 'w') as f:
        json.dump({
            'results': [asdict(topic) for topic in results.topics],
            'statistics': stats
        }, f, indent=2)
    
    # Print summary statistics
    print("\nSearch Method Comparison Results:")
    print("=" * 40)
    for method, rate in stats.items():
        if method.endswith('_success_rate'):
            print(f"{method}: {rate:.2f}%")
        elif method.endswith('_avg_time'):
            print(f"{method}: {rate:.2f} seconds")
        elif method.endswith('_unique_successes'):
            print(f"{method}: {int(rate)} unique successes")
    
    # Assert minimum success rate for best method or combination
    best_rate = max(rate for method, rate in stats.items() if method.endswith('_success_rate'))
    assert best_rate >= 50, f"Best success rate {best_rate:.2f}% is below minimum threshold of 50%"

def test_pamam_search():
    # Initialize the agent
    config = LLMModelConfig(
        provider=ProviderType.OPENAI,
        model_name="gpt-3.5-turbo",
        api_key=None
    )
    llm_service = LLMService(config)
    agent = PubChemAgent(llm_service)
    
    # Test different variations of the name
    variations = [
        "PAMAM dendrimer",
        "PAMAM-G0",
        "PAMAM G0",
        "Starburst PAMAM dendrimer",
        "Generation 0 PAMAM dendrimer"
    ]
    
    for query in variations:
        print(f"\nTesting search with: {query}")
        print("=" * 40)
        
        # Try direct search
        try:
            results = agent._search_pubchem_direct(query)
            if results:
                print(f"Direct search found: CID {results[0].cid}")
            else:
                print("Direct search found no results")
        except Exception as e:
            print(f"Direct search error: {str(e)}")
            
        # Try REST search
        try:
            results = agent._search_pubchem_rest(query)
            if results:
                print(f"REST search found: CID {results[0].cid}")
            else:
                print("REST search found no results")
        except Exception as e:
            print(f"REST search error: {str(e)}")
            
        # Try normalized search
        try:
            normalized = agent._normalize_query(query)
            print(f"Normalized variations: {normalized}")
            results = agent._search_with_fallbacks(query)
            if results:
                print(f"Fallback search found: CID {results[0].cid}")
            else:
                print("Fallback search found no results")
        except Exception as e:
            print(f"Fallback search error: {str(e)}")
            
def test_pamam_comparison():
    # Initialize the agent
    config = LLMModelConfig(
        provider=ProviderType.OPENAI,
        model_name="gpt-3.5-turbo",
        api_key=None
    )
    llm_service = LLMService(config)
    agent = PubChemAgent(llm_service)
    
    # Terms to compare
    terms = ["PAMAM", "PAMAM dendrimer"]
    
    for term in terms:
        print(f"\nResults for: {term}")
        print("=" * 50)
        
        # Test each search method
        methods = {
            'direct_rest_exact': agent._search_pubchem_direct,
            'direct_rest_fuzzy': agent._search_pubchem_rest,
            'pubchempy_name': agent._search_pubchem_direct,
            'pubchempy_formula': agent._search_pubchem_rest,
            'sourcetable': agent._search_pubchem_direct
        }
        
        for method_name, method_func in methods.items():
            print(f"\nMethod: {method_name}")
            print("-" * 30)
            
            start_time = time.time()
            try:
                results = method_func(term)
                processing_time = time.time() - start_time
                
                if results:
                    print(f"Success! Found {len(results)} results")
                    for i, result in enumerate(results[:3], 1):  # Show first 3 results
                        print(f"Result {i}:")
                        print(f"  CID: {result.cid}")
                        print(f"  Name: {result.name if hasattr(result, 'name') else 'N/A'}")
                        print(f"  Formula: {result.molecular_formula if hasattr(result, 'molecular_formula') else 'N/A'}")
                else:
                    print("No results found")
                print(f"Processing time: {processing_time:.2f}s")
                
            except Exception as e:
                print(f"Error: {str(e)}")
                print(f"Processing time: {time.time() - start_time:.2f}s")

if __name__ == "__main__":
    test_pamam_comparison() 