from agent_management.agents.pubchem_agent import PubChemAgent
from agent_management.llm_service import LLMService, LLMModelConfig, ProviderType
import json
import pytest
from typing import Dict, List, Tuple
import logging
import datetime
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
import concurrent.futures
import time
import sys

# Configure retries for requests
retry_strategy = Retry(
    total=3,
    backoff_factor=1,
    status_forcelist=[429, 500, 502, 503, 504],
)
adapter = HTTPAdapter(max_retries=retry_strategy)
http = requests.Session()
http.mount("https://", adapter)
http.mount("http://", adapter)

CHEMISTRY_TOPICS = [
    "Teach me about cryptands and 3D host-guest complexation",
    "Teach me about norbornane and its bridging structure",
    "Teach me about carboranes and their polyhedral cage structures",
    "Teach me about the DNA double helix and base-pair stacking",
    "Teach me about metal-metal quadruple bonds in dimolybdenum complexes",
    "Teach me about rotaxanes and their mechanically interlocked architecture",
    "Teach me about PAMAMdendrimers and their hyperbranched growth patterns",
    "Teach me about polyoxometalates and their metal-oxygen clusters",
    "Teach me about alkali-doped fullerenes and superconductivity",
    
    "Teach me about helicenes and their helical chirality",
    "Teach me about metallophthalocyanines and their planar macrocycles",
    "Teach me about organosilanes and the silicon hypervalency debate",
    "Teach me about zintl clusters and their electron-rich frameworks",
    "Teach me about peptidic β-sheets and α-helices in proteins",
    "Teach me about bridging metal-carbonyl ligands in cluster compounds",
    "Teach me about the Jahn-Teller distortion in octahedral Cu(II) complexes",
    "Teach me about cuneane and its strained cage system",
    "Teach me about tetraphenylporphyrin and its planar macrocycle",
    "Teach me about the Kekulé structure of benzene",
    "Teach me about the Woodward-Hoffmann rules for conrotatory and disrotatory reactions",
    "Teach me about circulene and its unique ring structure",
    "Teach me about diiron nonacarbonyl Fe2(CO)9 and its bridging CO groups",
    "Teach me about cryptophanes and their host-guest chemistry",
    "Teach me about double helical sulfur (S∞ chains) and polysulfur rings",
    "Teach me about the Schrock carbene complexes and metal-ligand multiple bonds",
    "Teach me about hydrogen-bonded molecular knots and trefoil structures",
    "Teach me about buckminsterfullerene (C60)",
    "Teach me about catenanes and how their rings interlock",
];


"""
this test should intake the chemistry topics and output the success results at the top of the json (percentage) and then the results for each topic below
PubChemAgent.get_molecule_sdfs() is like this: def get_molecule_sdfs(self, user_input: str, skip_llm: bool = False) -> PubChemSearchResult:
        1. Uses LLM to interpret the user input into a molecule name or ID.
        2. Queries PubChem for up to 5 matches.
        3. Fetches each match's SDF.
        4. Returns a structured result with the data.
        if skip_llm:
            # If skip_llm is True, directly use the input as the molecule name
            molecule_name = user_input
        else:
            # Step 1: Interpret the user's input via LLM
            molecule_name = self.interpret_user_query(user_input)
            secondary_molecule_name = self.interpret_user_query(user_input+" but not "+molecule_name + "so try to find a similar molecule that would be relevant to the user's query")
            if molecule_name == "N/A":
                # If the LLM can't parse or guess a name, return an empty structure
                return PubChemSearchResult(
                    query=user_input,
                    interpreted_query=molecule_name,
                    results=[]
                )

"""

def process_topic_with_timeout(agent: PubChemAgent, topic: str, timeout: int = 30) -> Dict:
    """
    Process a single topic with timeout.
    
    Args:
        agent: Initialized PubChemAgent
        topic: Topic to process
        timeout: Timeout in seconds
        
    Returns:
        Dictionary containing the result
    """
    try:
        with concurrent.futures.ThreadPoolExecutor() as executor:
            future = executor.submit(agent.get_molecule_sdfs, topic)
            result = future.result(timeout=timeout)
            
            success = len(result.results) > 0
            return {
                "topic": topic,
                "success": success,
                "interpreted_molecule": result.interpreted_query,
                "num_molecules_found": len(result.results),
                "molecules": [
                    {
                        "name": mol.name,
                        "cid": mol.cid,
                        "formula": mol.molecular_formula
                    } for mol in result.results
                ] if result.results else [],
                "error": None
            }
    except concurrent.futures.TimeoutError:
        return {
            "topic": topic,
            "success": False,
            "interpreted_molecule": None,
            "num_molecules_found": 0,
            "molecules": [],
            "error": "Request timed out"
        }
    except Exception as e:
        return {
            "topic": topic,
            "success": False,
            "interpreted_molecule": None,
            "num_molecules_found": 0,
            "molecules": [],
            "error": str(e)
        }

def run_chemistry_topics(topics: List[str], agent: PubChemAgent) -> Tuple[float, List[Dict], List[Dict]]:
    """
    Run tests on a list of chemistry topics and return success rate and results.
    
    Args:
        topics: List of topics to test
        agent: Initialized PubChemAgent
        
    Returns:
        Tuple of (success_rate, successful_results, failed_results)
    """
    logger = logging.getLogger(__name__)
    successful_results = []
    failed_results = []
    successful_queries = 0
    
    # Process each topic with timeout
    for topic in topics:
        logger.info(f"Processing topic: {topic}")
        result = process_topic_with_timeout(agent, topic)
        
        if result["success"]:
            successful_queries += 1
            successful_results.append(result)
        else:
            failed_results.append(result)
            logger.warning(f"Failed to process topic '{topic}': {result.get('error')}")
            
        # Add a small delay between requests to avoid rate limiting
        time.sleep(1)
    
    # Calculate success rate
    total_topics = len(topics)
    success_rate = (successful_queries / total_topics) * 100
    
    return success_rate, successful_results, failed_results

def test_chemistry_topics():
    """Test PubChemAgent's ability to handle various chemistry topics."""
    
    # Setup logging
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)
    
    # Initialize LLM service
    llm_config = LLMModelConfig(
        provider=ProviderType.ANTHROPIC,
        model_name="claude-3-5-sonnet-latest",
        max_tokens=1000
    )
    llm_service = LLMService(llm_config)
    
    # Initialize PubChemAgent
    agent = PubChemAgent(llm_service)
    
    # Run tests on all topics
    success_rate, successful_results, failed_results = run_chemistry_topics(CHEMISTRY_TOPICS, agent)
    
    # Prepare final report
    report = {
        "summary": {
            "success_rate": success_rate,
            "total_topics": len(CHEMISTRY_TOPICS),
            "successful_queries": len(successful_results),
            "failed_queries": len(failed_results),
            "timestamp": datetime.datetime.now().isoformat()
        },
        "successful_results": successful_results,
        "failed_results": failed_results
    }
    
    # Save results to file
    with open("pubchem_test_results.json", "w") as f:
        json.dump(report, f, indent=2)
        
    # Save failed queries to a separate file for focused testing
    with open("pubchem_failed_queries.json", "w") as f:
        json.dump({
            "timestamp": datetime.datetime.now().isoformat(),
            "failed_queries": [result["topic"] for result in failed_results]
        }, f, indent=2)
        
    logger.info(f"Success rate: {success_rate:.2f}%")
    logger.info(f"Failed queries saved to pubchem_failed_queries.json")
    
    # Assert minimum success rate
    assert success_rate >= 95, f"Success rate {success_rate:.2f}% is below target of 95%"
    
    return report

def test_failed_queries():
    """Test PubChemAgent specifically on previously failed queries."""
    
    try:
        with open("pubchem_failed_queries.json", "r") as f:
            failed_data = json.load(f)
            failed_queries = failed_data["failed_queries"]
    except FileNotFoundError:
        pytest.skip("No failed queries file found. Run test_chemistry_topics first.")
        return
        
    if not failed_queries:
        pytest.skip("No failed queries to test.")
        return
        
    # Initialize services
    llm_config = LLMModelConfig(
        provider=ProviderType.ANTHROPIC,
        model_name="claude-3-5-sonnet-latest",
        max_tokens=1000
    )
    llm_service = LLMService(llm_config)
    agent = PubChemAgent(llm_service)
    
    # Run tests on failed queries
    success_rate, successful_results, still_failing = run_chemistry_topics(failed_queries, agent)
    
    # Save focused test results
    report = {
        "summary": {
            "original_failures": len(failed_queries),
            "now_successful": len(successful_results),
            "still_failing": len(still_failing),
            "improvement_rate": success_rate,
            "timestamp": datetime.datetime.now().isoformat()
        },
        "now_successful": successful_results,
        "still_failing": still_failing
    }
    
    with open("pubchem_failed_queries_retest.json", "w") as f:
        json.dump(report, f, indent=2)
        
    # Log results
    logger.info(f"Retested {len(failed_queries)} previously failed queries")
    logger.info(f"Now successful: {len(successful_results)}")
    logger.info(f"Still failing: {len(still_failing)}")
    logger.info(f"Improvement rate: {success_rate:.2f}%")
    
    return report

def test_single_query():
    """Test PubChemAgent with a single query to verify basic functionality and timing."""
    
    # Setup detailed logging
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    logger = logging.getLogger(__name__)
    
    # Initialize services
    llm_config = LLMModelConfig(
        provider=ProviderType.ANTHROPIC,
        model_name="claude-3-5-sonnet-latest",
        max_tokens=1000
    )
    llm_service = LLMService(llm_config)
    agent = PubChemAgent(llm_service)
    
    # Test query
    test_query = "Teach me about water molecules and their 3D structure"
    logger.info(f"Starting test with query: {test_query}")
    
    # Record start time
    start_time = time.time()
    
    try:
        # Process with different timeouts to see what works
        timeouts = [30, 60, 120]  # Try 30s, 60s, 120s
        success = False
        
        for timeout in timeouts:
            try:
                logger.info(f"Attempting query with {timeout}s timeout...")
                result = process_topic_with_timeout(agent, test_query, timeout)
                
                # Log the result
                logger.info(f"Query completed in {time.time() - start_time:.2f} seconds")
                logger.info(f"Success: {result['success']}")
                logger.info(f"Interpreted molecule: {result['interpreted_molecule']}")
                logger.info(f"Number of molecules found: {result['num_molecules_found']}")
                
                if result['molecules']:
                    logger.info("Found molecules:")
                    for mol in result['molecules']:
                        logger.info(f"  - Name: {mol['name']}")
                        logger.info(f"    CID: {mol['cid']}")
                        logger.info(f"    Formula: {mol['formula']}")
                
                if result['success']:
                    success = True
                    break
                
            except concurrent.futures.TimeoutError:
                logger.warning(f"Timeout occurred after {timeout} seconds")
                continue
            except Exception as e:
                logger.error(f"Error during processing: {str(e)}")
                break
        
        # Save detailed results
        detailed_result = {
            "query": test_query,
            "total_time": time.time() - start_time,
            "success": success,
            "result": result if success else None,
            "timestamp": datetime.datetime.now().isoformat()
        }
        
        with open("single_query_test_result.json", "w") as f:
            json.dump(detailed_result, f, indent=2)
            
        logger.info("Test results saved to single_query_test_result.json")
        
        # Assert test passed
        assert success, "Single query test failed"
        
    except Exception as e:
        logger.error(f"Test failed with error: {str(e)}")
        raise

def test_all_topics_with_timeout():
    """Test all chemistry topics using the new timeout-based search method with detailed logging."""
    
    # Setup detailed logging
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    logger = logging.getLogger(__name__)
    
    # Initialize services
    llm_config = LLMModelConfig(
        provider=ProviderType.ANTHROPIC,
        model_name="claude-3-5-sonnet-latest",
        max_tokens=1000
    )
    llm_service = LLMService(llm_config)
    agent = PubChemAgent(llm_service)
    
    # Initialize results
    all_results = []
    successful_topics = 0
    total_topics = len(CHEMISTRY_TOPICS)
    start_time_all = time.time()
    
    logger.info(f"Starting test of {total_topics} topics")
    
    try:
        for index, topic in enumerate(CHEMISTRY_TOPICS, 1):
            topic_start_time = time.time()
            
            # Show progress
            progress = (index / total_topics) * 100
            sys.stdout.write(f"\rProcessing topic {index}/{total_topics} ({progress:.1f}%)")
            sys.stdout.flush()
            
            logger.info(f"\nProcessing topic {index}/{total_topics}: {topic}")
            
            # Process topic with timeout
            try:
                result = process_topic_with_timeout(agent, topic, timeout=30)  # Use 30s timeout
                
                # Add timing information
                result["processing_time"] = time.time() - topic_start_time
                
                # Log detailed results
                logger.info(f"Topic completed in {result['processing_time']:.2f} seconds")
                logger.info(f"Success: {result['success']}")
                logger.info(f"Interpreted molecule: {result['interpreted_molecule']}")
                
                if result["success"]:
                    successful_topics += 1
                    logger.info("Found molecules:")
                    for mol in result["molecules"]:
                        logger.info(f"  - Name: {mol['name']}")
                        logger.info(f"    CID: {mol['cid']}")
                        logger.info(f"    Formula: {mol['formula']}")
                else:
                    logger.warning(f"Failed to process topic: {result.get('error')}")
                
                all_results.append(result)
                
            except Exception as e:
                logger.error(f"Error processing topic: {str(e)}")
                all_results.append({
                    "topic": topic,
                    "success": False,
                    "error": str(e),
                    "processing_time": time.time() - topic_start_time
                })
            
            # Add delay between requests
            time.sleep(1)
        
        # Calculate final statistics
        total_time = time.time() - start_time_all
        success_rate = (successful_topics / total_topics) * 100
        
        # Prepare comprehensive report
        report = {
            "summary": {
                "total_topics": total_topics,
                "successful_topics": successful_topics,
                "failed_topics": total_topics - successful_topics,
                "success_rate": success_rate,
                "total_time": total_time,
                "average_time_per_topic": total_time / total_topics,
                "timestamp": datetime.datetime.now().isoformat()
            },
            "results": all_results
        }
        
        # Save detailed report
        with open("timeout_test_results.json", "w") as f:
            json.dump(report, f, indent=2)
        
        # Save failed topics for focused testing
        failed_topics = [r["topic"] for r in all_results if not r["success"]]
        with open("timeout_failed_topics.json", "w") as f:
            json.dump({
                "timestamp": datetime.datetime.now().isoformat(),
                "failed_topics": failed_topics
            }, f, indent=2)
        
        # Log final summary
        logger.info("\n" + "="*50)
        logger.info("Test Summary:")
        logger.info(f"Total topics processed: {total_topics}")
        logger.info(f"Successful topics: {successful_topics}")
        logger.info(f"Failed topics: {total_topics - successful_topics}")
        logger.info(f"Success rate: {success_rate:.2f}%")
        logger.info(f"Total time: {total_time:.2f} seconds")
        logger.info(f"Average time per topic: {(total_time / total_topics):.2f} seconds")
        logger.info("="*50)
        
        # Assert minimum success rate
        assert success_rate >= 95, f"Success rate {success_rate:.2f}% is below target of 95%"
        
        return report
        
    except Exception as e:
        logger.error(f"Test failed with error: {str(e)}")
        raise

def test_multiple_queries():
    """Test PubChemAgent with the first five topics to verify functionality and timing."""
    
    # Setup detailed logging
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    logger = logging.getLogger(__name__)
    
    # Initialize services
    llm_config = LLMModelConfig(
        provider=ProviderType.ANTHROPIC,
        model_name="claude-3-5-sonnet-latest",
        max_tokens=1000
    )
    llm_service = LLMService(llm_config)
    agent = PubChemAgent(llm_service)
    
    # Take first 5 topics
    test_topics = CHEMISTRY_TOPICS[:5]
    results = []
    successful_topics = 0
    start_time_all = time.time()
    
    logger.info(f"Starting test with {len(test_topics)} topics")
    
    try:
        for index, topic in enumerate(test_topics, 1):
            topic_start_time = time.time()
            logger.info(f"\nTesting topic {index}/{len(test_topics)}: {topic}")
            
            # Process with 30s timeout (since this worked in the single query test)
            try:
                result = process_topic_with_timeout(agent, topic, timeout=30)
                
                # Add timing information
                result["processing_time"] = time.time() - topic_start_time
                
                # Log the result
                logger.info(f"Query completed in {result['processing_time']:.2f} seconds")
                logger.info(f"Success: {result['success']}")
                logger.info(f"Interpreted molecule: {result['interpreted_molecule']}")
                logger.info(f"Number of molecules found: {result['num_molecules_found']}")
                
                if result['success']:
                    successful_topics += 1
                    logger.info("Found molecules:")
                    for mol in result['molecules']:
                        logger.info(f"  - Name: {mol['name']}")
                        logger.info(f"    CID: {mol['cid']}")
                        logger.info(f"    Formula: {mol['formula']}")
                
                results.append(result)
                
            except Exception as e:
                logger.error(f"Error processing topic: {str(e)}")
                results.append({
                    "topic": topic,
                    "success": False,
                    "error": str(e),
                    "processing_time": time.time() - topic_start_time
                })
            
            # Add delay between requests
            time.sleep(1)
        
        # Calculate final statistics
        total_time = time.time() - start_time_all
        success_rate = (successful_topics / len(test_topics)) * 100
        
        # Save detailed results
        detailed_result = {
            "summary": {
                "total_topics": len(test_topics),
                "successful_topics": successful_topics,
                "failed_topics": len(test_topics) - successful_topics,
                "success_rate": success_rate,
                "total_time": total_time,
                "average_time_per_topic": total_time / len(test_topics),
                "timestamp": datetime.datetime.now().isoformat()
            },
            "results": results
        }
        
        with open("multiple_queries_test_result.json", "w") as f:
            json.dump(detailed_result, f, indent=2)
            
        logger.info("\n" + "="*50)
        logger.info("Test Summary:")
        logger.info(f"Total topics processed: {len(test_topics)}")
        logger.info(f"Successful topics: {successful_topics}")
        logger.info(f"Failed topics: {len(test_topics) - successful_topics}")
        logger.info(f"Success rate: {success_rate:.2f}%")
        logger.info(f"Total time: {total_time:.2f} seconds")
        logger.info(f"Average time per topic: {(total_time / len(test_topics)):.2f} seconds")
        logger.info("="*50)
        
        return detailed_result
        
    except Exception as e:
        logger.error(f"Test failed with error: {str(e)}")
        raise